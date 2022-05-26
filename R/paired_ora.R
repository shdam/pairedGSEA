#' Paired Over-Representation Analysis
#' 
#' paired_ora uses \code{\link[fgsea:fora]{fora()}} to run the over-representation analysis.
#'   First the aggregated pvalues are adjusted using the Benjamini & Hochberg method.
#'   The analysis is run on all significant genes found by DESeq2 and DEXSeq individually.
#'   I.e., two runs of fora() are executed and subsequently joined into a single object.
#' 
#' @param paired_de_result The output of \code{\link[pairedGSEA:paired_de]{paired_de()}}
#' @param gene_sets List of gene sets to analyse
#' @param cutoff (Default: 0.05) Adjusted p-value cutoff for selecting significant genes
#' @param min_size (Default: 25) Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param experiment_title Title of your experiment. Your results will be stored in paste0("results/", experiment_title, "_fora.RDS").
#' @param quiet (Default: FALSE) Whether to print messages
#' @family paired
#' @export 
paired_ora <- function(paired_de_result,
                       gene_sets,
                       cutoff = 0.05,
                       min_size = 25,
                       experiment_title = NULL,
                       quiet = FALSE){
  
  # Check column names are as expected
  check_colname(paired_de_result, "pvalue_deseq", "paired_de_result")
  check_colname(paired_de_result, "pvalue_dexseq", "paired_de_result")
  check_colname(paired_de_result, "gene", "paired_de_result")
  
  # Significant genes
  if(!quiet)  message("Identifying differentially expressed genes")
  genes_deseq <- paired_de_result %>% 
    dplyr::filter(!is.na(pvalue_deseq) & !is.na(gene),
                  padj_deseq < cutoff) %>%
    dplyr::arrange(padj_deseq)
  
  genes_dexseq <- paired_de_result %>% 
    dplyr::filter(!is.na(pvalue_dexseq) & !is.na(gene),
                  padj_dexseq < cutoff) %>%
    dplyr::arrange(padj_dexseq)
  
  # fora
  if(!quiet) message("Running over-representation analysis")
  universe <- unique(paired_de_result$gene)
  ## ORA on DESeq2 results
  ora_deseq <- fgsea::fora(gene_sets, genes = genes_deseq$gene, 
                           universe = universe, minSize = min_size) %>% 
    dplyr::rename(size_geneset = size) %>% 
    dplyr::mutate(size_genes = nrow(genes_deseq),
                  size_universe = length(universe),
                  relative_risk = (overlap / size_geneset) / (size_genes / size_universe),
                  enrichment_score = log2(relative_risk + 0.06)) 
  
  ## ORA on DEXSeq results
  ora_dexseq <- fgsea::fora(gene_sets, genes = genes_dexseq$gene,
                            universe = universe, minSize = min_size) %>% 
    dplyr::rename(size_geneset = size) %>% 
    dplyr::mutate(size_genes = nrow(genes_dexseq),
                  size_universe = length(universe),
                  relative_risk = (overlap / size_geneset) / (size_genes / size_universe),
                  enrichment_score = log2(relative_risk + 0.06))
  
  if(!quiet) message("Joining result")
  ora_joined <- ora_deseq %>% 
    dplyr::full_join(ora_dexseq,
                     by = c("pathway", "size_geneset", "size_universe"),
                     suffix = c("_deseq", "_dexseq")) %>% 
    dplyr::mutate(enrichment_score_shift = relative_risk_dexseq - relative_risk_deseq,
                  enrichment_score_shift = log2(abs(enrichment_score_shift)) * sign(enrichment_score_shift),
                  experiment = experiment_title) %>% 
    dplyr::arrange(enrichment_score_shift)
  
  if(!is.null(experiment_title)){
    if(!quiet) message("Storing fora results")
    store_result(ora_joined, paste0(experiment_title, "_ora.RDS"), "ORA on both DESeq2 and DEXSeq results", quiet = quiet)
  }
  
  
  return(ora_joined)
}

