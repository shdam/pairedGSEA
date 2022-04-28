#' Run paired DESeq2 and DEXSeq analyses
#' @param metadata A metadata file or data frame object
#' @param group_col The metadata column specifying the what group each sample is associated with
#' @param sample_col The column in the metadata that specifies the sample IDs (should correspond to column names in tx_count)
#' @param comparison The comparison to use for this particular experiment. Format example: "1v2". In that example "1" would define the baseline and "2" would define to case samples.
#' @param tx_count The transcripts count matrix of an RNA-seq analysis
#' @param experiment_title Title of your experiment. Your results will be stored in paste0("results/", experiment_title, "_pairedGSEA.RDS").
#' @param run_sva (Default: TRUE) A logical stating whether SVA should be run.
#' @param prefilter (Default: 10) The prefilter threshold, where rowSums lower than the prefilter threshold will be removed from the count matrix. Set to 0 or FALSE to prevent prefiltering
#' @param fit_type (Default: "local") Either "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of dispersions to the mean intensity.
#' @param store_results (Default: TRUE) A logical indicating if results should be stored in the folder "results/".
#' @param quiet (Default: FALSE) Whether to print messages
#' @param parallel (Default: FALSE) If FALSE, no parallelization. If TRUE, parallel execution using BiocParallel, see next argument BPPARAM.
#' @param BPPARAM (Default: \code{BiocParallel::bpparam()}) An optional parameter object passed internally to bplapply when parallel = TRUE. If not specified, the parameters last registered with register will be used.
#' 
#' @export
paired_gsea <- function(tx_count,
                        metadata,
                        group_col,
                        sample_col,
                        comparison,
                        experiment_title = "MyLittleExperiment",
                        store_results = TRUE,
                        run_sva = TRUE,
                        prefilter = 10,
                        fit_type = "local",
                        quiet = FALSE,
                        parallel = FALSE,
                        BPPARAM = BiocParallel::bpparam()){
  
  # Check parallel
  if(parallel) check_missing_package("BiocParallel", "Bioc")
  
  # Ensure correct format for comparison
  comparison <- check_comparison(comparison)
  
  if(!quiet) message("Running ", experiment_title)
  
  # Load metadata
  if(!quiet) message("Preparing metadata")
  metadata <- prepare_metadata(metadata, group_col, comparison)
  
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
  dds <- convert_matrix_to_dds(tx_count, metadata, group_col)
  
  if(run_sva){
    dds <- run_sva(dds, group_col, quiet = quiet)
  }
  


  ### Run DESeq2
  deseq_results <- run_deseq(
    dds,
    group_col = group_col,
    comparison = comparison,
    fit_type = fit_type,
    experiment_title = experiment_title,
    store_results = store_results,
    quiet = quiet,
    parallel = parallel,
    BPPARAM = BPPARAM
    )
  
  # Store results
  if(store_results) store_result(deseq_results, paste0(experiment_title, "_deseq2res.RDS"), "DESeq2 results", quiet = quiet)
  
  
  ### Run DEXSeq
  dexseq_results <- run_dexseq(
    dds,
    group_col = group_col,
    comparison = comparison,
    experiment_title = experiment_title,
    store_results = store_results,
    quiet = quiet,
    parallel = parallel,
    BPPARAM = BPPARAM
    )
  # Store results
  if(store_results) store_result(dexseq_results, paste0(experiment_title, "_dexseqres.RDS"), "DEXSeq results", quiet = quiet)
  
  if(!quiet) message(experiment_title, " is analysed.")
  
  return(list("deseq" = deseq_results, "dexseq" = dexseq_results))
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
#' @inheritParams paired_gsea
#' @keywords internal
prepare_metadata <- function(metadata, group_col, comparison){
  
  if(typeof(metadata[1]) == "character" & length(metadata) == 1){
    if(stringr::str_ends(metadata, ".xlsx")) metadata <- readxl::read_excel(metadata)
    else if(stringr::str_ends(metadata, ".csv")) metadata <- readr::read_csv(metadata)
  } else if("data.frame" %!in% class(metadata)){stop("Please provide path to a metadata file or a data.frame object.")}
  
  if(group_col %!in% colnames(metadata)) stop("Could not find column ", group_col, " in metadata.")
  # Ensure comparison is on the right format
  comparison <- check_comparison(comparison)
  # Remove irrelevant groups
  metadata <- metadata[metadata[[group_col]] %in% comparison, ]
  
  # Check metadata content
  if(nrow(metadata) == 0) stop("The values in", group_col, "does not correspond to the comparison values.\nPlease ensure what comes before and after the 'v' in comparison are what is found in the", group_col, "column.")
  
  # Add comparison levels to metadata
  metadata[[group_col]] <- factor(metadata[[group_col]], levels = comparison)
  
  return(metadata)
}




#' Run SVA on DESeqDataSet
#' 
#' This internal function runs a surrogate variable analysis on the count matrix in the DESeqDataSet. 
#'   The found surrogate variables will then be added to the metadata and the design formula in the DESeqDataSet object to be used in the DGE and DTU analyses.
#'  
#' @inheritParams paired_gsea
#' @param dds A DESeqDataSet. See \code{?DESeq2::DESeqDataSet} for more information about the object type.
#' @keywords internal
run_sva <- function(dds, group_col, quiet = FALSE){
  
  
  if(!quiet) message("Running SVA")
  # Normalize counts with DESeq2 for SVA
  normalized_counts <- DESeq2::normTransform(dds) %>% 
    SummarizedExperiment::assay()
  
  # Extract metadata
  metadata <- SummarizedExperiment::colData(dds)
  # Define model matrix 
  mod1 <- stats::model.matrix(~metadata[[group_col]])
  mod0 <- cbind(mod1[, 1])
  
  # Run SVA
  svseq <- sva::sva(normalized_counts, mod1, mod0)
  message("Found ", svseq$n.sv, " surrogate variables")
  # Store surrogate variables and rename for ease of reference
  svs <- tibble::as_tibble(svseq$sv, .name_repair = "minimal")
  colnames(svs) <- paste0("sv", 1:svseq$n.sv)
  
  if(!quiet) message("Redefining DESeq design formula\n")
  # Add svs to dds colData
  SummarizedExperiment::colData(dds) <- cbind(metadata, svs)
  # Redefine design formula to include svs
  DESeq2::design(dds) <- stats::formula(paste0("~", group_col, "+", stringr::str_c(colnames(svs), collapse = "+")))
  
  
  return(dds)
}


#' Run DEXSeq analysis
#' 
#' This internal function runs a differential transcript usage analysis using DEXSeq (See \code{?DEXSeq::DEXSeq} for more detailed information).
#'   Here, the surrogate variables found by \code{\link{run_sva}}, if any, will be added to the DEXSeqDataSet before running the analysis.
#' 
#' @inheritParams paired_gsea
#' @inheritParams run_sva
#' @keywords internal
run_dexseq <- function(dds,
                       group_col,
                       comparison,
                       experiment_title = "MyLittleExperiment",
                       store_results = FALSE,
                       quiet = FALSE,
                       parallel = FALSE,
                       BPPARAM = BiocParallel::bpparam()){
  
  if(parallel) check_missing_package("BiocParallel", "Bioc")
  
  if(!quiet) message("Initiating DEXSeq")
  # Ensure correct format for comparison
  comparison <- check_comparison(comparison)
  
  # Extract group and feature from rownames of DESeq2 object
  group_feat <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  
  # Extract the found surrogate variables
  svs <- as.character(DESeq2::design(dds))[2] %>% 
    stringr::str_split(" \\+ ", n = 2, simplify = TRUE) %>% 
    .[2] %>% 
    stringr::str_split(" \\+ ", simplify = TRUE)
  
  # Add surrogate variables to DEXSeq design formula
  if(svs == ""){
    new_design <- "~sample + exon + condition:exon"
    reduced_design <- "~sample + exon"
  } else{
    new_design <- paste0("~sample + exon + condition:exon + ", stringr::str_c(svs, ":exon", collapse = " + "))
    reduced_design <- paste0("~sample + exon + ", stringr::str_c(svs, ":exon", collapse = " + "))
  }
  design_formula <- stats::formula(
    new_design
  )
  reduced_formula <- stats::formula(
    reduced_design
  )
  
  # Define sample data based on DESeq2 object
  sample_data <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% 
    dplyr::select(dplyr::all_of(group_col), dplyr::starts_with("sv")) %>% 
    dplyr::rename(condition = dplyr::all_of(group_col)) %>% 
    dplyr::mutate(condition = dplyr::case_when(condition == comparison[1] ~ "B",      # baseline 
                                               condition == comparison[2] ~ "C") %>%  # condition
                    factor(levels = c("C", "B")))
  
  
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
    dplyr::rename(log2FC_dexseq = log2fold_B_C)
  
  return(dexseq_results)
}







#' Run DESeq2 analysis
#' @inheritParams paired_gsea
#' @inheritParams run_sva
#' @keywords internal
run_deseq <- function(dds,
                      group_col,
                      comparison,
                      fit_type = "local",
                      experiment_title = "MyLittleExperiment",
                      store_results = FALSE,
                      quiet = FALSE,
                      parallel = FALSE,
                      BPPARAM = BiocParallel::bpparam()){
  
  # Register parallel
  if(parallel) check_missing_package("BiocParallel", "Bioc")
  
  if(!quiet) message("Running DESeq2")
  dds <- DESeq2::DESeq(dds, parallel = parallel, BPPARAM = BPPARAM,
                       fitType = fit_type, quiet = FALSE)
  
  # Store DESeqDataSet with DESeq2 analysis
  if(store_results) {
    store_result(dds, paste0(experiment_title, "_dds.RDS"), "DESeqDataSet", quiet = quiet)
  }
  
  # Ensure correct format for comparison
  comparison <- check_comparison(comparison)
  
  if(!quiet) message("Extracting results")
  deseq_results <- DESeq2::results(dds,
                                   contrast = c(group_col, comparison),
                                   parallel = parallel,
                                   BPPARAM = BPPARAM)
  
  
  # Convert result to tibble
  deseq_results <- deseq_results %>% 
    tibble::as_tibble(rownames = "gene_tx") %>% 
    tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
    dplyr::rename(log2FC_deseq = log2FoldChange)
  
  return(deseq_results)
}
