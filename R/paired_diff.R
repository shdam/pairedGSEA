#' Run paired DESeq2 and DEXSeq analyses
#' 
#' With paired_diff you can run a paired differential gene expression and
#' splicing analysis. The function expects a counts matrix or a
#' SummarizedExperiment or DESeqDataSet object as input.
#' A preliminary prefiltering step is performed to remove genes with a summed 
#' count lower than the provided threshold. Likewise, genes with counts in 
#' only one sample are removed. This step is mostly to speed up differential 
#' analyses, as DESeq2 will do a stricter filtering.
#' Surrogate Variable Analysis is recommended to allow the analyses to take
#' batch effects etc. into account.
#' After the two differential analyses, the transcript-level p-values will be
#' aggregated to gene-level to allow subsequent Gene-Set Enrichment Analysis.
#' Transcript-level results can be extracted by setting
#' \code{store_results = TRUE}.
#' 
#' 
#' 
#' @inheritParams DESeq2::DESeq
#' @param object A data object of the types matrix, SummarizedExperiment,
#' or DESeqDataSet. If a matrix is used, please also provide metadata.
#' @param metadata (Default: NULL) A metadata file or data frame object
#' @param group_col The metadata column specifying the what group each
#' sample is associated with
#' @param sample_col The column in the metadata that specifies the sample IDs
#' (should correspond to column names in tx_count)
#' @param baseline Group value of baseline samples
#' @param case Group value of case samples
#' @param covariates Name of column(s) in the metadata that indicate(s)
#' covariates. E.g., c("gender", "tissue_type")
#' @param experiment_title Title of your experiment. Your results will be
#' stored in paste0("results/", experiment_title, "_pairedGSEA.RDS").
#' @param run_sva (Default: TRUE) A logical stating whether SVA should be run.
#' @param prefilter (Default: 10) The prefilter threshold, where rowSums lower
#' than the prefilter threshold will be removed from the count matrix.
#' Set to 0 or FALSE to prevent prefiltering
#' @param fit_type (Default: "local") Either "parametric", "local", "mean", or
#' "glmGamPoi" for the type of fitting of dispersions to the mean intensity.
#' @param store_results (Default: TRUE) A logical indicating if results should
#' be stored in the folder "results/".
#' @param quiet (Default: FALSE) Whether to print messages
#' @param expression_only (Default: FALSE) 
#' A logical that indicates whether to only
#' run DESeq2 analysis. Not generally recommended.
#' The setting was implemented to make the SVA impact analysis easier
#' @param custom_design (Default: FALSE) A logical or formula. Can be used to
#' apply a custom desing formula for the analysis. Generally not recommended, 
#' as pairedGSEA will make its own design formula
#' from the group and covariate columns
#' @param parallel (Default: FALSE) If FALSE, no parallelization.
#' If TRUE, parallel execution using BiocParallel, see next argument BPPARAM.
#' @param BPPARAM (Default: \code{BiocParallel::bpparam()}) An optional
#' parameter object passed internally to bplapply when parallel = TRUE.
#' If not specified, the parameters last registered with registerwill be used.
#' @param ... Additional parameters passed to
#' \code{\link[DESeq2:DESeq]{DESeq()}}
#' @family paired
#' @importFrom methods is
#' @import SummarizedExperiment
#' @import DESeq2
#' @import DEXSeq
#' @import BiocParallel
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate across all_of full_join
#' @importFrom stringr str_detect
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
#' @return A data.frame of aggregated pvalues
#' @usage 
#' paired_diff(
#'     object,
#'     group_col,
#'     sample_col,
#'     baseline,
#'     case,
#'     metadata = NULL,
#'     covariates = NULL,
#'     experiment_title = NULL,
#'     store_results = TRUE,
#'     run_sva = TRUE,
#'     prefilter = 10,
#'     test = "LRT",
#'     fit_type = "local",
#'     quiet = FALSE,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     expression_only = FALSE,
#'     custom_design = FALSE,
#'     ...
#'     )
#' @examples 
#' 
#' # Run analysis on included example data
#' data("example_se")
#' 
#' paired_diff(
#'     object = example_se,
#'     group_col = "group_nr",
#'     sample_col = "id",
#'     baseline = 1,
#'     case = 2,
#'     experiment_title = "Example",
#'     store_results = FALSE 
#' )
#' 
#' @export
paired_diff <- function(
        object,
        group_col,
        sample_col,
        baseline,
        case,
        metadata = NULL,
        covariates = NULL,
        experiment_title = NULL,
        store_results = TRUE,
        run_sva = TRUE,
        prefilter = 10,
        test = "LRT",
        fit_type = "local",
        quiet = FALSE,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        expression_only = FALSE,
        custom_design = FALSE,
        ...){

    ## Initial error checks
    stopifnot(
        "Please provide a valid data type: matrix, SummarizedExperiment,
        DESeqDataSet" = any(class(object) %in% c(
            "matrix", "SummarizedExperiment","DESeqDataSet")))
    
    stopifnot(
        "Cannot store results if experiment title haven't been given" =
            (store_results & !is.null(experiment_title)) | !store_results)
    
    stopifnot(
        "Please provide metadata with your count matrix" =
            (is(object, "matrix") & !is.null(metadata)) | !is(object, "matrix"))
    
    ## Checking column names
    stopifnot(
        "Covariate names must not contain spaces" =
            stringr::str_detect(covariates, " ", negate = TRUE))
    stopifnot(
        "group_col name must not contain spaces" =
            stringr::str_detect(group_col, " ", negate = TRUE))
    stopifnot(
        "sample_col name must not contain spaces" =
            stringr::str_detect(sample_col, " ", negate = TRUE))

    if(!quiet) message("Running ", experiment_title)

    ## Define design formula
    if(custom_design == TRUE & is(object, "DESeqDataSet")) {
        design <- DESeq2::design(object)
        } else if(is(custom_design, "formula")) {
            design <- custom_design
            } else {
                design <- formularise_vector(c(group_col, covariates))
                if(is(object, "DESeqDataSet")) warning(
                    "OBS: your design will be overwritten to: ",
                    as.character(design))
                }

    ## Convert se to dds
    if(is(object, "SummarizedExperiment")) {
        SummarizedExperiment::colData(object) <- SummarizedExperiment::colData(
            object) %>% 
            tibble::as_tibble() %>% 
            # Ensure columns are factors
            dplyr::mutate(dplyr::across(dplyr::all_of(
                c(group_col, covariates)), factor)) %>% 
            S4Vectors::DataFrame(row.names = colnames(object))
        object <- DESeq2::DESeqDataSet(object, design)
        }

    ## Load metadata
    if(!quiet) message("Preparing metadata")
    if(!(is(object, "matrix")) & is.null(metadata))
        metadata <- SummarizedExperiment::colData(object)

    metadata <- prepare_metadata(metadata, group_col, paste(c(baseline, case)))

    ## Subsample metadata to only include samples present in the count matrix
    metadata <- metadata[metadata[[sample_col]] %in% colnames(object), ]
    stopifnot(
        "Please ensure that the sample IDs in the metadata matches the
        column names of the count matrix." = nrow(metadata) > 0)
    
    ## Check for presence of undesired characters
    if(any(stringr::str_detect(metadata[[sample_col]], "[- ]")))
        message("OBS! Some or all sample names contain a '-' or ' ', ",
                "which will cause downstream methods to complain.")

    ## Check sample_col is in metadata
    stopifnot(
        "Sample column not in metadata" = sample_col %in% colnames(metadata))

    # Add metadata to DESeqDataSet
    if(is(object, "DESeqDataSet"))
        SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(metadata)

    ## Convert count matrix to DESeqDataSet
    if(is(object, "matrix")){
        ## Ensure rows in metadata matches columns in the count matrix
        object <- object[, metadata[[sample_col]]]
        
        ## Create DDS from count matrix
        if(!quiet) message("Converting count matrix to DESeqDataSet")
        object <- convert_matrix_to_dds(object, metadata, design)
        }

    stopifnot(
        "Please ensure the rownames have the format 'gene:transcript'" =
            (expression_only | (!expression_only & stringr::str_detect(
                rownames(object)[1], ":"))))

    # Rename object variable to dds
    dds <- object; rm(object)

    # Prefiltering
    if(prefilter) dds <- pre_filter(dds, prefilter)

    # Detect surrogate variables
    if(run_sva){
        dds <- run_sva(dds, quiet = quiet)
    }

    # Run DESeq2
    expression_results <- run_deseq(
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

    ## If only DGE is requested (used for SVA evaluations in the paper)
    if(expression_only){
        expression_aggregated <- aggregate_pvalue(
            expression_results, gene = "gene", weights = "baseMean",
            lfc = "lfc", type = "expression") %>% 
            dplyr::mutate(padj = stats::p.adjust(pvalue, "fdr"))

    if(store_results) store_result(
        expression_aggregated, paste0(
            experiment_title,
            "_aggregated_pvals.RDS"),
        "gene pvalue aggregation")
    return(expression_aggregated)
    }

    # Run DEXSeq
    splicing_results <- run_dexseq(
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


    # Aggregate p values
    if(!quiet) {
        if(!is(experiment_title, "NULL"))
            message(experiment_title, " is analysed.")
        message("Aggregating p values")
    }

    splicing_aggregated <- aggregate_pvalue(
        splicing_results, gene = "groupID", weights = "exonBaseMean",
        type = "splicing")
    expression_aggregated <- aggregate_pvalue(
        expression_results, gene = "gene", weights = "baseMean",
        type = "expression")

    aggregated_pvals <- dplyr::full_join(
        expression_aggregated,
        splicing_aggregated,
        by = "gene",
        suffix = c("_expression", "_splicing"))


    if(store_results) store_result(
        aggregated_pvals, paste0(experiment_title, "_aggregated_pvals.RDS"),
        "gene pvalue aggregation")

    if(!quiet) message("Done.")


    return(aggregated_pvals)
}


#' Prepare metadata
#' 
#' This internal function reads a filepath for a metadata file
#' (or a data.frame object), loads it, 
#' and filters out the rows that do not contain data from the
#' relevant patient group, as defined in the \code{comparison} argument.
#' 
#' @examples 
#' \dontrun{
#' metadata <- prepare_metadata(metadata = "path/to/metadata.xlsx",
#'     group_col = "group", comparison = "2v1")
#' }
#' 
#' @inheritParams paired_diff
#' @param baseline_case A character vector with baseline and case values
#' @importFrom stringr str_ends
#' @note Suggested: importFrom readxl read_excel importFrom readr read_csv
#' @keywords internal
#' @return A data.frame
prepare_metadata <- function(metadata, group_col, baseline_case){

    if(is(metadata, "character")){
        if(stringr::str_ends(metadata, ".xlsx")) {
            check_missing_package(package = "readxl")
            metadata <- readxl::read_excel(metadata)}
        else if(stringr::str_ends(metadata, ".csv")) {
            check_missing_package(package = "readr")
            metadata <- readr::read_csv(metadata)}
    } else if(!is(metadata, "data.frame") & !is(metadata, "DataFrame")){
            stop(
            "Please provide path to a metadata file or a data.frame /
            DataFrame object.")}

    if(group_col %!in% colnames(metadata))
        stop("Could not find column ", group_col, " in metadata.")
    
    # Remove irrelevant groups
    metadata <- metadata[
        as.character(metadata[[group_col]]) %in% baseline_case,]
    
    # Check both baseline and case is in metadata
    stopifnot(
        "Check for misspellings in baseline or case IDs." = 
            all(baseline_case %in% metadata[[group_col]])
        )

    # Check metadata content
    if(nrow(metadata) == 0) stop(
        "The values in", group_col,
        "does not correspond to the baseline and case values.")

    # Add comparison levels to metadata
    metadata[[group_col]] <- factor(
        metadata[[group_col]], levels = baseline_case)

    return(metadata)
}


#' Convert matrix to DESeqDataSet
#' 
#' This internal function converts a matrix and metadata object into a
#' DESeqDataSet, using the group_col as the basic design.
#' 
#' @param tx_count Count matrix to convert to dds
#' @inheritParams paired_diff
#' @keywords internal
#' @return A DESeqDataSet
convert_matrix_to_dds <- function(tx_count, metadata, design){

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = tx_count,
        colData = metadata,
        design = design)
    return(dds)
}

#' Run SVA on DESeqDataSet
#' 
#' This internal function runs a surrogate variable analysis on the
#' count matrix in the DESeqDataSet. 
#' The found surrogate variables will then be added to the metadata
#' and the design formula in the DESeqDataSet object to be used in
#' the DGE and DTU analyses.
#'  
#' @inheritParams paired_diff
#' @param dds A DESeqDataSet. See \code{?DESeq2::DESeqDataSet} 
#' for more information about the object type.
#' @importFrom stats model.matrix cor na.omit
#' @importFrom sva sva
#' @importFrom stringr str_split
#' @keywords internal
#' @return A DESeqDataSet
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
    if(!quiet) message("\nFound ", svseq$n.sv, " surrogate variables")

    if(svseq$n.sv == 0) return(dds)

    # Store surrogate variables and rename for ease of reference
    svs <- tibble::as_tibble(svseq$sv, .name_repair = "minimal")
    colnames(svs) <- paste0("sv", seq_len(svseq$n.sv))

    # Remove svs that confound with mod1
    cors <- as.matrix(stats::cor(svs, mod1[,2:ncol(mod1)]))
    cors[abs(cors) > 0.8] <- NA
    if(any(is.na(cors)) & !quiet) message(
        "Removing surrogate variables confounding with design model.")
    cors <- stats::na.omit(cors)
    svs <- svs[, rownames(cors)]

    if(!quiet) message("Redefining DESeq2 design formula\n")
    # Add svs to dds colData
    SummarizedExperiment::colData(dds) <- cbind(metadata, svs)
    # Redefine design formula to include svs
    DESeq2::design(dds) <- formularise_vector(c(
        as.character(design)[2] %>% 
            stringr::str_split(" \\+ ", simplify = TRUE), colnames(svs)))

    return(dds)
}


#' Run DESeq2 analysis
#' 
#' This internal function runs a differential gene expression analysis using
#' DESeq2 (See \code{?DESeq2::DESeq} for more detailed information).
#' Here, the surrogate variables found by \code{\link{run_sva}},
#' if any, will be added to the DESeqDataSet before running the analysis.
#'   
#' @inheritParams paired_diff
#' @inheritParams run_sva
#' @param ... Additional parameters passed to
#' \code{\link[DESeq2:DESeq]{DESeq()}}
#' @importFrom tidyr separate
#' @importFrom dplyr rename
#' @keywords internal
#' @return A data.frame
#' @usage
#' run_deseq(
#'     dds,
#'     group_col,
#'     baseline,
#'     case,
#'     test = "LRT",
#'     fit_type = "local",
#'     experiment_title = "Experiment Title",
#'     store_results = FALSE,
#'     quiet = FALSE,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     ...
#'     )
#' 
run_deseq <- function(
        dds,
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


    dds <- DESeq2::DESeq(
        dds,
        reduced = reduced,
        test = test,
        parallel = parallel, 
        BPPARAM = BPPARAM,
        fitType = fit_type,
        quiet = quiet,
        ...)

    # Store DESeqDataSet with DESeq2 analysis
    if(store_results) {
        store_result(
            dds, paste0(experiment_title, "_dds.RDS"),
            "DESeqDataSet", quiet = quiet)
    }


    if(!quiet) message("Extracting results")
    expression_results <- DESeq2::results(
        dds,
        contrast = c(group_col, paste(c(case, baseline))),
        parallel = parallel,
        BPPARAM = BPPARAM)

    # Store result extraction
    if(store_results) {
        if(store_results) store_result(
            expression_results,
            paste0(experiment_title, "_expression_results.RDS"),
            "differential expression results", quiet = quiet)
    }

    # Convert result to tibble
    expression_results <- expression_results %>% 
        tibble::as_tibble(rownames = "gene_tx") %>% 
        tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
        dplyr::rename(lfc = log2FoldChange)

    return(expression_results)
}


#' Run DEXSeq analysis
#' 
#' This internal function runs a differential transcript usage analysis using
#' DEXSeq (See \code{?DEXSeq::DEXSeq} for more detailed information).
#' Here, the surrogate variables found by \code{\link{run_sva}}, if any,
#' will be added to the DEXSeqDataSet before running the analysis.
#' 
#' @inheritParams paired_diff
#' @inheritParams run_sva
#' @importFrom dplyr select all_of rename mutate case_when
#' @keywords internal
#' @return A data.frame
#' @usage
#' run_dexseq(
#'     dds,
#'     group_col,
#'     baseline,
#'     case,
#'     experiment_title = "Experiment Title",
#'     store_results = FALSE,
#'     quiet = FALSE,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam()
#'     )
#' 
run_dexseq <- function(
        dds,
        group_col,
        baseline,
        case,
        experiment_title = "Experiment Title",
        store_results = FALSE,
        quiet = FALSE,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam()){

    if(!quiet) message("Initiating DEXSeq")

    stopifnot(
        "Please ensure the rownames have the format 'gene:transcript'" =
            stringr::str_detect(rownames(dds)[1], ":"))

    # Extract group and feature from rownames of DESeq2 object
    group_feat <- rownames(dds) %>% 
        stringr::str_split(":", simplify = TRUE)

    # Extract the found surrogate variables and  covariates
    svs_covariates <- DESeq2::design(dds) %>% 
        reduce_formula(formularise = FALSE)

    # Add surrogate variables and covariates to DEXSeq design formula
    if(svs_covariates[[1]] == "1"){
        svs_covariates <- NULL
        design_formula <- formularise_vector(
            c("sample", "exon", "condition:exon"))
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
        as.data.frame(row.names = SummarizedExperiment::colnames(dds)) %>% 
        dplyr::select(dplyr::all_of(c(group_col, svs_covariates))) %>% 
        dplyr::rename(condition = dplyr::all_of(group_col)) %>% 
        dplyr::mutate(condition = dplyr::case_when(
            condition == baseline ~ "B", 
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
        store_result(
            dxd, paste0(experiment_title, "_dxd.RDS"), "DEXSeqDataSet",
            quiet = quiet)
    }

    ### Run DEXSeq
    if(!quiet) message("\nRunning DEXSeq -- This might take a while")
    if(!parallel) BiocParallel::register(BiocParallel::SerialParam())
    splicing_results <- DEXSeq::DEXSeq(
        dxd,
        reducedModel = reduced_formula,
        BPPARAM = BPPARAM,
        quiet = quiet)

    # Store result extraction
    if(store_results) store_result(
        splicing_results, paste0(experiment_title, "_splicing_results.RDS"),
        "differnetial splicing results", quiet = quiet)

    # Rename LFC column for consistency and human-readability
    splicing_results <- splicing_results %>% 
        tibble::as_tibble() %>% 
        dplyr::rename(lfc = log2fold_C_B)

    return(splicing_results)
}




#' p-value aggregation to gene level
#' 
#' This internal function aggregates p-values to gene level using
#' Lancaster aggregation.
#' 
#' @param df A data.frame-type object with pvalues to aggregate
#' @param gene (Default: "gene") The column with gene names
#' @param p (Default: "pvalue") The column with p-values
#' @param lfc (Default: "lfc") The column with log2-foldchanges
#' @param type (Default: "expression") 
#' Either "expression" for aggregation of
#' DESeq2 results or "splicing" for aggregation of DEXSeq results
#' @param weights (Default: "baseMean") The column to use
#' for weighting the aggregation
#' @importFrom dplyr filter rename mutate group_by summarise case_when
#' @importFrom aggregation lancaster
#' @importFrom stats weighted.mean
#' @keywords internal
#' @usage 
#' aggregate_pvalue(
#'     df,
#'     gene = "gene",
#'     p = "pvalue",
#'     weights = "baseMean",
#'     lfc = "lfc",
#'     type = "expression"
#'     )
#' @return A data.frame
aggregate_pvalue <- function (
        df,
        gene = "gene",
        p = "pvalue",
        weights = "baseMean",
        lfc = "lfc",
        type = "expression"){
    stopifnot("df is not a data.frame-type object" = is(df, "data.frame"))
    stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))

    type <- tolower(type)

    res <- df %>% 
        dplyr::filter(!is.na(padj)) %>% 
        dplyr::rename(
            pvalue = .data[[p]],
            gene = .data[[gene]],
            lfc = .data[[lfc]],
            weights = .data[[weights]]) %>% 
        # Prevent warning from Lancaster
        dplyr::mutate(pvalue = dplyr::case_when(
            pvalue < 10e-320 ~ 10e-320,
            TRUE ~ pvalue)) %>% 
        dplyr::group_by(gene)
    # p value aggregation
    if(type == "splicing"){
        res <- res %>% 
            dplyr::summarise(
            lfc = lfc[which(pvalue == min(pvalue))][[1]],
            pvalue = aggregation::lancaster(pvalue, weights)
            )
        } else if(type == "expression"){
        res <- res %>% 
            dplyr::summarise(
                lfc = stats::weighted.mean(lfc, weights),
                pvalue = aggregation::lancaster(pvalue, weights)
                )
        }
    # Adjust p values and remove zeros to prevent downstream issues
    res <- res %>% 
        dplyr::mutate(
            padj = stats::p.adjust(pvalue, "fdr"),
            pvalue = dplyr::case_when(
                pvalue < 10e-320 ~ 10e-320,
                TRUE ~ pvalue)) 
    return(res)
}