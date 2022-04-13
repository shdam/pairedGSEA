#' Run paired DESeq2 and DEXSeq analyses
#' @param metadata A metadata file or data frame object
#' @param group_col The metadata column specifying the what group each sample is associated with
#' @param sample_col The column in the metadata that specifies the sample IDs (should correspond to column names in tx_count)
#' @param comparison The comparison to use for this particular experiment. Format example: "1v2"
#' @param tx_count The transcripts count matrix of an RNA-seq analysis
#' @param experiment_title Title of your experiment. Your results will be stored in paste0("results/", experiment_title, "_pairedGSEA.RDS").
#' @param run_sva (Default: TRUE) A logical stating whether SVA should be run.
#' @param prefilter (Default: 10) The prefilter threshold, where rowSums lower than the prefilter threshold will be removed from the count matrix. Set to 0 or FALSE to prevent prefiltering
#' @param fit_type (Default: "local") Either "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of dispersions to the mean intensity.
#' @param dds_out (Default: FALSE) Put a file string to store DESeqDataSet after DESeq analysis, before results are extracted
#' @param dxd_out (Default: FALSE) Put a file string to store DEXSeqDataSet after DEXSeq analysis, before results are extracted
#' @param quiet (Default: FALSE) Whether to print messages
#' @param parallel (Default: FALSE) If FALSE, no parallelization. If TRUE, parallel execution using BiocParallel, see next argument BPPARAM.
#' @param BPPARAM (Default: \code{BiocParallel::bpparam()}) An optional parameter object passed internally to bplapply when parallel = TRUE. If not specified, the parameters last registered with register will be used.
#' 
#' @importFrom magrittr %>%
#' @export
paired_gsea <- function(tx_count,
                        metadata,
                        group_col,
                        sample_col,
                        comparison,
                        experiment_title = "MyLittleExperiment",
                        run_sva = TRUE,
                        prefilter = 10,
                        fit_type = "local",
                        dds_out = FALSE,
                        dxd_out = FALSE,
                        quiet = FALSE,
                        parallel = FALSE,
                        BPPARAM = BiocParallel::bpparam()){
  
  # Check parallel
  if(parallel) pairedGSEA:::check_missing_package("BiocParallel", "Bioc")
  
  # Ensure correct format for comparison
  comparison <- pairedGSEA:::check_comparison(comparison)
  
  if(!quiet) message("Running ", experiment_title)
  
  # Load metadata
  if(!quiet) message("Preparing metadata")
  metadata <- pairedGSEA:::prepare_metadata(metadata, group_col, comparison)
  
  # Check sample_col is in metadata
  stopifnot("Sample column not in metadata" = sample_col %in% colnames(m))
  
  # Subsample metadata to only include samples present in the count matrix
  metadata <- metadata[metadata[[sample_col]] %in% colnames(tx_count), ]
  stopifnot("Please ensure that the sample IDs in the metadata matches the column names of the count matrix." = nrow(metadata) > 0)
  # Ensure rows in metadata matches columns in the count matrix
  tx_count <- tx_count[, metadata[[sample_col]]]
  
  ### Prefiltering
  if(prefilter) tx_count <- pairedGSEA:::pre_filter(tx_count, prefilter)
  
  if(!quiet) message("Converting count matrix to DESeqDataSet")
  # Create DDS from count matrix
  dds <- pairedGSEA:::convert_matrix_to_dds(tx_count, metadata, group_col)
  
  if(run_sva){
    dds <- pairedGSEA:::run_sva(dds, group_col, quiet = quiet)
  }
  

  ### Run DESeq2
  deseq_results <- pairedGSEA:::run_deseq(
    dds,
    group_col = group_col,
    comparison = comparison,
    fit_type = fit_type,
    dds_out = dds_out,
    quiet = quiet,
    parallel = parallel,
    BPPARAM = BiocParallel::bpparam()
    )
  
  # Store results
  pairedGSEA:::store_result(deseq_results, paste0(experiment_title, "_deseq2res.RDS"), "DESeq2 results", quiet = quiet)
  rm(deseq_results)
  
  ### Run DEXSeq
  dexseq_results <- pairedGSEA:::run_dexseq(
    dds,
    group_col = group_col,
    comparison = comparison,
    dxd_out = dxd_out,
    quiet = quiet,
    parallel = parallel,
    BPPARAM = BiocParallel::bpparam()
    )
  # Store results
  pairedGSEA:::store_result(dexseq_results, paste0(experiment_title, "_dexseqres.RDS"), "DEXSeq results", quiet = quiet)
  
  if(!quiet) message(experiment_title, " is analysed.")
}
