
#' Run DESeq2
#' @inheritParams paired_gsea
run_deseq <- function(dds,
                      group_col,
                      comparison,
                      fit_type = "local",
                      dds_out = FALSE,
                      quiet = FALSE,
                      parallel = FALSE,
                      BPPARAM = BiocParallel::bpparam()){
  
  # Register parallel
  if(parallel) pairedGSEA:::check_missing_package("BiocParallel", "Bioc")
  
  if(!quiet) message("Running DESeq2")
  dds <- DESeq2::DESeq(dds, parallel = parallel, BPPARAM = BPPARAM,
                       fitType = fit_type, quiet = FALSE)
  
  # Store DESeqDataSet with DESeq2 analysis
  if(typeof(dds_out) == "character") {
    pairedGSEA:::store_result(dds, dds_out, "DESeqDataSet", quiet = quiet)
  }
  
  # Ensure correct format for comparison
  comparison <- pairedGSEA:::check_comparison(comparison)

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
