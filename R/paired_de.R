#' Run paired DESeq2 and DEXSeq analyses
#' 
#' @inheritParams DESeq2::DESeq
#' @param metadata A metadata file or data frame object
#' @param group_col The metadata column specifying the what group each sample is associated with
#' @param sample_col The column in the metadata that specifies the sample IDs (should correspond to column names in tx_count)
#' @param baseline Group value of baseline samples
#' @param case Group value of case samples
#' @param covariates Name of column(s) in the metadata that indicate(s) covariates. E.g., c("gender", "tissue_type")
#' @param tx_count The transcripts count matrix of an RNA-seq analysis
#' @param experiment_title Title of your experiment. Your results will be stored in paste0("results/", experiment_title, "_pairedGSEA.RDS").
#' @param run_sva (Default: TRUE) A logical stating whether SVA should be run.
#' @param prefilter (Default: 10) The prefilter threshold, where rowSums lower than the prefilter threshold will be removed from the count matrix. Set to 0 or FALSE to prevent prefiltering
#' @param fit_type (Default: "local") Either "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of dispersions to the mean intensity.
#' @param store_results (Default: TRUE) A logical indicating if results should be stored in the folder "results/".
#' @param quiet (Default: FALSE) Whether to print messages
#' @param parallel (Default: FALSE) If FALSE, no parallelization. If TRUE, parallel execution using BiocParallel, see next argument BPPARAM.
#' @param BPPARAM (Default: \code{BiocParallel::bpparam()}) An optional parameter object passed internally to bplapply when parallel = TRUE. If not specified, the parameters last registered with register will be used.
#' @param ... Additional parameters passed to \code{\link[DESeq2:DESeq]{DESeq()}}
#' @family paired
#' @export
paired_de <- function(tx_count,
                      metadata,
                      group_col,
                      sample_col,
                      baseline,
                      case,
                      covariates = NULL,
                      experiment_title = "Experiment Title",
                      store_results = TRUE,
                      run_sva = TRUE,
                      prefilter = 10,
                      test = "LRT",
                      fit_type = "local",
                      quiet = FALSE,
                      parallel = FALSE,
                      BPPARAM = BiocParallel::bpparam(),
                      ...){
  
  # Checking column names
  stopifnot("Covariate names must not contain spaces" = stringr::str_detect(covariates, " ", negate = TRUE))
  stopifnot("group_col name must not contain spaces" = stringr::str_detect(group_col, " ", negate = TRUE))
  
  if(!quiet) message("Running ", experiment_title)
  
  # Load metadata
  if(!quiet) message("Preparing metadata")
  metadata <- prepare_metadata(metadata, group_col, paste(c(baseline, case)))
  
  # Check sample_col is in metadata
  stopifnot("Sample column not in metadata" = sample_col %in% colnames(metadata))
  
  # Subsample metadata to only include samples present in the count matrix
  metadata <- metadata[metadata[[sample_col]] %in% colnames(tx_count), ]
  stopifnot("Please ensure that the sample IDs in the metadata matches the column names of the count matrix." = nrow(metadata) > 0)
  
  # Ensure rows in metadata matches columns in the count matrix
  tx_count <- tx_count[, metadata[[sample_col]]]
  
  ### Prefiltering
  if(prefilter) tx_count <- pre_filter(tx_count, prefilter)
  
  if(!quiet) message("Converting count matrix to DESeqDataSet")
  # Create DDS from count matrix
  design <- formularise_vector(c(group_col, covariates))
  dds <- convert_matrix_to_dds(tx_count, metadata, design)
  
  if(run_sva){
    dds <- run_sva(dds, quiet = quiet)
  }
  


  ### Run DESeq2
  deseq_results <- run_deseq(
    dds,
    group_col = group_col,
    baseline = baseline,
    case = case,
    fit_type = fit_type,
    test = test,
    experiment_title = experiment_title,
    store_results = store_results,
    quiet = quiet,
    parallel = parallel,
    BPPARAM = BPPARAM,
    ...
    )
  
  # Store results
  if(store_results) store_result(deseq_results, paste0(experiment_title, "_deseqres.RDS"), "DESeq2 results", quiet = quiet)
  
  
  ### Run DEXSeq
  dexseq_results <- run_dexseq(
    dds,
    group_col = group_col,
    baseline = baseline,
    case = case,
    experiment_title = experiment_title,
    store_results = store_results,
    quiet = quiet,
    parallel = parallel,
    BPPARAM = BPPARAM
    )
  
  # Store results
  if(store_results) store_result(dexseq_results, paste0(experiment_title, "_dexseqres.RDS"), "DEXSeq results", quiet = quiet)
  
  if(!quiet) message(experiment_title, " is analysed.")
  
  if(!quiet) message("Aggregating p values")
  
  dexseq_aggregated <- aggregate_pvalue(dexseq_results, gene = "groupID", weights = "exonBaseMean", lfc = "log2FC_dexseq", type = "dexseq")
  deseq_aggregated <- aggregate_pvalue(deseq_results, gene = "gene", weights = "baseMean", lfc = "log2FC_deseq", type = "deseq")
  
  aggregated_pvals <- dplyr::full_join(deseq_aggregated,
                                       dexseq_aggregated,
                                       by = "gene",
                                       suffix = c("_deseq", "_dexseq")) %>% 
    dplyr::mutate(padj_deseq = stats::p.adjust(pvalue_deseq, "fdr"),
                  padj_dexseq = stats::p.adjust(pvalue_dexseq, "fdr"))
  
  
  
  if(store_results) store_result(aggregated_pvals, paste0(experiment_title, "_aggregated_pvals.RDS"), "gene pvalue aggregation")
  
  
  return(aggregated_pvals)
}


#' Prepare metadata
#' 
#' This internal function reads a filepath for a metadata file (or a data.frame object), loads it, 
#'   and filters out the rows that do not contain data from the relevant patient group, as defined in the \code{comparison} argument.
#' 
#' @examples 
#' \dontrun{
#' metadata <- prepare_metadata(metadata = "path/to/metadata.xlsx",
#'   group_col = "group", comparison = "2v1")
#' }
#' 
#' @inheritParams paired_de
#' @param baseline_case A character vector with baseline and case values
#' @keywords internal
prepare_metadata <- function(metadata, group_col, baseline_case){
  
  if(typeof(metadata[1]) == "character" & length(metadata) == 1){
    if(stringr::str_ends(metadata, ".xlsx")) metadata <- readxl::read_excel(metadata)
    else if(stringr::str_ends(metadata, ".csv")) metadata <- readr::read_csv(metadata)
  } else if("data.frame" %!in% class(metadata)){stop("Please provide path to a metadata file or a data.frame object.")}
  
  if(group_col %!in% colnames(metadata)) stop("Could not find column ", group_col, " in metadata.")

  # Remove irrelevant groups
  metadata <- metadata[as.character(metadata[[group_col]]) %in% baseline_case, ]
  
  # Check metadata content
  if(nrow(metadata) == 0) stop("The values in", group_col, "does not correspond to the baseline and case values.")
  
  # Add comparison levels to metadata
  metadata[[group_col]] <- factor(metadata[[group_col]], levels = baseline_case)
  
  return(metadata)
}


#' Convert matrix to DESeqDataSet
#' 
#' This internal function converts a matrix and metadata object into a DESeqDataSet, using the group_col as the basic design.
#' 
#' @inheritParams paired_de
#' @keywords internal
convert_matrix_to_dds <- function(tx_count, metadata, design){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = tx_count,
                                        colData = metadata,
                                        design = design)
  return(dds)
}

#' Run SVA on DESeqDataSet
#' 
#' This internal function runs a surrogate variable analysis on the count matrix in the DESeqDataSet. 
#'   The found surrogate variables will then be added to the metadata and the design formula in the DESeqDataSet object to be used in the DGE and DTU analyses.
#'  
#' @inheritParams paired_de
#' @param dds A DESeqDataSet. See \code{?DESeq2::DESeqDataSet} for more information about the object type.
#' @keywords internal
run_sva <- function(dds, quiet = FALSE){
  
  
  if(!quiet) message("Running SVA")
  # Normalize counts with DESeq2 for SVA
  normalized_counts <- DESeq2::normTransform(dds) %>% 
    SummarizedExperiment::assay()
  
  # Extract metadata
  metadata <- SummarizedExperiment::colData(dds)
  
  # Get dds design
  design <- DESeq2::design(dds)
  reduced <- design %>% 
    reduce_formula()
    
  
  # Define model matrix 
  mod1 <- stats::model.matrix(design, data = metadata)
  mod0 <- stats::model.matrix(reduced, data = metadata)
  
  # Run SVA
  svseq <- sva::sva(normalized_counts, mod = mod1, mod0 = mod0)
  if(!quiet) message("Found ", svseq$n.sv, " surrogate variables")
  
  if(svseq$n.sv == 0) return(dds)
  
  # Store surrogate variables and rename for ease of reference
  svs <- tibble::as_tibble(svseq$sv, .name_repair = "minimal")
  colnames(svs) <- paste0("sv", 1:svseq$n.sv)
  
  if(!quiet) message("Redefining DESeq design formula\n")
  # Add svs to dds colData
  SummarizedExperiment::colData(dds) <- cbind(metadata, svs)
  # Redefine design formula to include svs
  DESeq2::design(dds) <- formularise_vector(c(as.character(design)[2], colnames(svs)))
  
  return(dds)
}


#' Run DESeq2 analysis
#' 
#' This internal function runs a differential gene expression analysis using DESeq2 (See \code{?DESeq2::DESeq} for more detailed information).
#'   Here, the surrogate variables found by \code{\link{run_sva}}, if any, will be added to the DESeqDataSet before running the analysis.
#'   
#' @inheritParams paired_de
#' @inheritParams run_sva
#' @param ... Additional parameters passed to \code{\link[DESeq2:DESeq]{DESeq()}}
#' @keywords internal
run_deseq <- function(dds,
                      group_col,
                      baseline,
                      case,
                      test = "LRT",
                      fit_type = "local",
                      experiment_title = "Experiment Title",
                      store_results = FALSE,
                      quiet = FALSE,
                      parallel = FALSE,
                      BPPARAM = BiocParallel::bpparam(),
                      ...){
  
  if(!quiet) message("Running DESeq2")
  
  # Reduce design formula to surrogate variables and covariates
  reduced <- DESeq2::design(dds) %>% 
    reduce_formula()

  
  dds <- DESeq2::DESeq(dds,
                       reduced = reduced,
                       test = test,
                       parallel = parallel, 
                       BPPARAM = BPPARAM,
                       fitType = fit_type,
                       quiet = quiet,
                       ...)
  
  # Store DESeqDataSet with DESeq2 analysis
  if(store_results) {
    store_result(dds, paste0(experiment_title, "_dds.RDS"), "DESeqDataSet", quiet = quiet)
  }
  
  
  if(!quiet) message("Extracting results")
  deseq_results <- DESeq2::results(dds,
                                   contrast = c(group_col, paste(c(case, baseline))),
                                   parallel = parallel,
                                   BPPARAM = BPPARAM)
  
  
  # Convert result to tibble
  deseq_results <- deseq_results %>% 
    tibble::as_tibble(rownames = "gene_tx") %>% 
    tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
    dplyr::rename(log2FC_deseq = log2FoldChange)
  
  return(deseq_results)
}


#' Run DEXSeq analysis
#' 
#' This internal function runs a differential transcript usage analysis using DEXSeq (See \code{?DEXSeq::DEXSeq} for more detailed information).
#'   Here, the surrogate variables found by \code{\link{run_sva}}, if any, will be added to the DEXSeqDataSet before running the analysis.
#' 
#' @inheritParams paired_de
#' @inheritParams run_sva
#' @keywords internal
run_dexseq <- function(dds,
                       group_col,
                       baseline,
                       case,
                       experiment_title = "Experiment Title",
                       store_results = FALSE,
                       quiet = FALSE,
                       parallel = FALSE,
                       BPPARAM = BiocParallel::bpparam()){
  
  if(!quiet) message("Initiating DEXSeq")
  
  # Extract group and feature from rownames of DESeq2 object
  group_feat <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  
  # Extract the found surrogate variables and  covariates
  svs_covariates <- DESeq2::design(dds) %>% 
    reduce_formula(formularise = FALSE)
  
  # Add surrogate variables and covariates to DEXSeq design formula
  if(svs_covariates[[1]] == "1"){
    design_formula <- formularise_vector(c("sample", "exon", "condition:exon"))
    reduced_formula <- formularise_vector(c("sample", "exon"))
  } else{
    design_formula <- c("sample", "exon", "condition:exon") %>% 
      c(paste0(svs_covariates, ":exon")) %>% 
      formularise_vector()
    reduced_formula <- c("sample", "exon") %>% 
      c(paste0(svs_covariates, ":exon")) %>% 
      formularise_vector()
  }
  
  # Define sample data based on DESeq2 object
  sample_data <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% 
    dplyr::select(dplyr::all_of(c(group_col, svs_covariates))) %>% 
    dplyr::rename(condition = dplyr::all_of(group_col)) %>% 
    dplyr::mutate(condition = dplyr::case_when(condition == baseline ~ "B", 
                                               condition == case ~ "C") %>% 
                    factor(levels = c("B", "C")))
  
  
  # Convert to DEXSeqDataSet
  if(!quiet) message("Creating DEXSeqDataSet")
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = DESeq2::counts(dds),
    sampleData = sample_data,
    design = design_formula,
    groupID = group_feat[, 1],
    featureID = group_feat[, 2]
  )
  
  # Store DEXSeqDataSet before DEXSeq analysis
  if(store_results) {
    store_result(dxd, paste0(experiment_title, "_dxd.RDS"), "DEXSeqDataSet", quiet = quiet)
  }
  
  ### Run DEXSeq
  if(!quiet) message("Running DEXSeq")
  if(!parallel) BiocParallel::register(BiocParallel::SerialParam())
  dexseq_results <- DEXSeq::DEXSeq(dxd,
                                   reducedModel = reduced_formula,
                                   BPPARAM = BPPARAM,
                                   quiet = quiet)
  
  
  # Rename LFC column for consistency and human-readability
  dexseq_results <- dexseq_results %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(log2FC_dexseq = log2fold_C_B)
  
  return(dexseq_results)
}




#' p-value aggregation to gene level
#' 
#' This internal function aggregates p-values to gene level using Lancaster agreggation..
#' 
#' @param df A data.frame-type object with pvalues to aggregate
#' @param gene (Default: "gene") The column with gene names
#' @param p (Default: "pvalue") The column with p-values
#' @param lfc (Default: "log2FC") The column with log2-foldchanges
#' @param type (Default: "deseq") Either "deseq" for aggregation of DESeq2 results or "dexseq" for aggregation of DEXSeq results
#' @param weights (Default: "baseMean") The column to use for weighting the aggregation
#' @keywords internal
aggregate_pvalue <- function (df,
                              gene = "gene",
                              p = "pvalue",
                              weights = "baseMean",
                              lfc = "log2FC",
                              type = "deseq"){
  stopifnot("df is not a data.frame-type object" = "data.frame" %in% class(df))
  stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))
  
  type <- tolower(type)
  
  res <- df %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::rename(pvalue = .data[[p]],
                  gene = .data[[gene]],
                  lfc = .data[[lfc]],
                  weights = .data[[weights]]) %>% 
    # Prevent warning from Lancaster
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) %>% 
    dplyr::group_by(gene) %>% 
    # p value aggregation
    purrr::when(type == "dexseq" ~ 
                  dplyr::summarise(.,
                                   lfc = lfc[which(pvalue == min(pvalue))][[1]],
                                   pvalue = aggregation::lancaster(pvalue, weights)),
                type == "deseq" ~ 
                  dplyr::summarise(.,
                                   lfc = weighted.mean(lfc, weights),
                                   pvalue = aggregation::lancaster(pvalue, weights))) %>% 
    # Remove zeros to prevent downstream issues
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) 
  return(res)
}