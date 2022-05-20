#' Paired Over-Representation Analysis
#' 
#' paired_ora uses \code{\link[fgsea:fora]{fora()}} to run the over-representation analysis.
#'   First the aggregated pvalues are adjusted using the Benjamini & Hochberg method.
#'   The analysis is run on all significant genes found by DESeq2 and DEXSeq individually.
#'   I.e., two runs of fora() are executed and subsequently joined into a single object.
#' 
#' @param paired_de_result The output of \code{\link[pairedGSEA:paired_de]{paired_de()}}
#' @param gene_sets List of gene sets to analyse
#' @param min_size (Default: 25) Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param experiment_title Title of your experiment. Your results will be stored in paste0("results/", experiment_title, "_fora.RDS").
#' @param quiet (Default: FALSE) Whether to print messages
#' @family paired
#' @export 
paired_ora <- function(paired_de_result,
                       gene_sets,
                       min_size = 25,
                       experiment_title = NULL,
                       quiet = FALSE){
  
  # Check column names are as expected
  check_colname(paired_de_result, "pvalue_deseq", "paired_de_result")
  check_colname(paired_de_result, "pvalue_dexseq", "paired_de_result")
  check_colname(paired_de_result, "ensembl_gene", "paired_de_result")
  
  # Significant genes
  if(!quiet)  message("Identifying differentially expressed genes")
  genes_deseq <- paired_de_result %>% 
    dplyr::filter(!is.na(pvalue_deseq) & !is.na(ensembl_gene)) %>%
    dplyr::mutate(padj = stats::p.adjust(pvalue_deseq, "fdr")) %>% 
    dplyr::filter(padj < 0.05) %>% 
    dplyr::arrange(padj)
  
  genes_dexseq <- paired_de_result %>% 
    dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>%
    dplyr::mutate(padj = stats::p.adjust(pvalue_dexseq, "fdr")) %>% 
    dplyr::filter(padj < 0.05) %>% 
    dplyr::arrange(padj)
  
  # fora
  if(!quiet) message("Running over-representation analysis")
  universe <- unique(paired_de_result$ensembl_gene)
  ## ORA on DESeq2 results
  fora_deseq <- fgsea::fora(gene_sets, genes = genes_deseq$ensembl_gene, 
                            universe = universe, minSize = min_size) %>% 
    dplyr::rename(size_geneset = size) %>% 
    dplyr::mutate(size_genes = nrow(genes_deseq),
                  size_universe = length(universe),
                  odds_ratio = (overlap / size_genes) / (size_geneset / size_universe),
                  enrichment_score = log2(odds_ratio)) 
  
  ## ORA on DEXSeq results
  fora_dexseq <- fgsea::fora(gene_sets, genes = genes_dexseq$ensembl_gene,
                             universe = universe, minSize = min_size) %>% 
    dplyr::rename(size_geneset = size) %>% 
    dplyr::mutate(size_genes = nrow(genes_dexseq),
                  size_universe = length(universe),
                  odds_ratio = (overlap / (size_geneset-overlap)) / (size_genes / (size_universe-size_genes)),
                  enrichment_score = log2(odds_ratio))
  
  if(!quiet) message("Joining result")
  fora_joined <- fora_deseq %>% 
    dplyr::full_join(fora_dexseq, by = c("gene_sets", "size_geneset", "size_universe"), suffix = c("_deseq", "_dexseq"))
  fora_joined$experiment_title <- experiment_title
  
  if(!is.null(experiment_title)){
    if(!quiet) message("Storing fora results")
    # pairedGSEA:::store_result(fora_deseq, paste0(experiment_title, "_fora_deseq.RDS"), "fora on DESeq2 results", quiet = quiet)
    # pairedGSEA:::store_result(fora_dexseq, paste0(experiment_title, "_fora_dexseq.RDS"), "fora on DEXSeq results", quiet = quiet)
    store_result(fora_joined, paste0(experiment_title, "_fora.RDS"), "fora on both DESeq2 and DEXSeq results", quiet = quiet)
  }
  
  
  return(fora_joined)
}


