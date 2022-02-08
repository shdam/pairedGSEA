
#' Run DESeq2 and DEXSeq analyses
#' 
#' @importFrom tidyr pivot_wider
#' @export
runExperiment <- function(row, archs4db = NULL, txCount = NULL, groupCol = "group_nr", tpm = TRUE, prefilter = 10, parallel = TRUE){
  
  if(typeof(row) == "character"){ # Convert apply-made row to tibble
    row <- tibble::as_tibble(row, rownames = "names") %>% 
      tidyr::pivot_wider(values_from = value, names_from = names)
  }
  
  
  message("Running on ", row$study)
  
  ### Load metadata
  md_file <- row$filename
  dataname <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove("csv") 
  
  ### Define tpm file
  if(tpm) tpm <- stringr::str_replace(archs4db, "counts", "tpm")
  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)`
  experimentTitle <- row$`comparison_title (empty_if_not_okay)`
  
  
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  
  
  ### Prepare for DE
  dds <- prepDE(md = md_file,
                gtf = gtf,
                archs4db = archs4db,
                txCount = txCount,
                groupCol = groupCol,
                comparison = comparison,
                prefilter = prefilter)
  

  
  ### Run DESeq2
  res <- runDESeq2(dds,
                   groupCol = groupCol,
                   comparison = comparison,
                   samples = dds$id,
                   tpm = tpm,
                   gtf = gtf,
                   parallel = parallel,
                   fitType = "local",
                   BPPARAM = BiocParallel::bpparam())#, dds_out = "deseq2_1_GSE154968.RDS")
  
  # Store results
  check_make_dir("results")
  saveRDS(res, paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
  rm(res)
  
  ### Run DEXSeq
  dxr <- runDEXSeq(dds, groupCol, comparison)
  # Store results
  saveRDS(dxr, paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
  
  message(row$study, " is done.")
  
}

#' Analyse experiments
#' 
#' @importFrom fgsea fgseaMultilevel
#' @export
analyseExperiment <- function(row){
  
  deseq2file <- paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS")
  if(!file.exists(deseq2file)) stop(paste0("File:", deseq2file, " does not exists.\\n","Please run experiment before analysing. See ?runExperiment"))
  
  if(typeof(row) == "character"){ # Convert apply-made row to tibble
    row <- tibble::as_tibble(row, rownames = "names") %>% 
      tidyr::pivot_wider(values_from = value, names_from = names)
  }
  
  
  ### Load metadata
  md_file <- row$filename
  dataname <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove("csv") 
  
  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)`
  experimentTitle <- row$`comparison_title (empty_if_not_okay)`
  
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  
  
  message("Analysing ", row$study, " ", experimentTitle)
  
  # Load results
  message("Loading results")
  
  res <- readRDS(paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
  dxr <- readRDS(paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
  
  # p value aggregation
  
  message("Aggregating p values")
  
  dxr_agg <- perGenePValue(dxr, gene = "groupID", weights = "exonBaseMean") %>% 
    dplyr::filter(pvalue < 0.05)
  res_agg <- perGenePValue(res, gene = "gene", weights = "baseMean", lfc = "log2FC") %>% 
    dplyr::filter(pvalue < 0.05)
  
  (comb <- dplyr::full_join(res_agg, dxr_agg, by = "ensembl_gene", suffix = c("_deseq2", "_dexseq")))
  
  
  message("Storing aggregation result")
  saveRDS(comb, paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
  
  
  message("Gene set enrichment analysis")
  
  ### Defining stats
  resStats <- comb %>% 
    dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_deseq2) * sign(lfc)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene)
  
  dxrStats <- comb %>% 
    dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_dexseq)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene) 
  
  ### Run fgsea
  fgseaRes <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = resStats,
                                     scoreType = "std",
                                     eps = 10e-320
  )
  saveRDS(fgseaRes, paste0("results/", dataname, "_fgseaRes_", experimentTitle, ".RDS"))
  # fgseaRes_sig <- fgseaRes %>% filter(padj<0.05)
  
  fgseaDxr <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = dxrStats,
                                     scoreType = "pos",
                                     eps = 10e-320
  )
  saveRDS(fgseaDxr, paste0("results/", dataname, "_fgseaDxr_", experimentTitle, ".RDS"))
  # fgseaDxr_sig <- fgseaDxr %>% filter(padj<0.05)
  
  message("fgsea results are stored in the results folder. Look for '*_fgseaRes_*' and '*_fgseaDxr_*'")
}

#' Per gene p value aggregation
#' 
#' @importFrom aggregation lancaster
#' @importFrom purrr when
#' @export
perGenePValue <- function (df,
                           gene = "gene",
                           p = "pvalue",
                           weights = "baseMean",
                           lfc = NULL){
  stopifnot(typeof(df) == "list")
  stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))
  
  res <- df %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::rename(pvalue = .data[[p]],
                  ensembl_gene = .data[[gene]]) %>% 
    # Prevent warning from Lancaster
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) %>% 
    dplyr::group_by(ensembl_gene) %>% 
    # p value aggregation
    purrr::when(is.null(lfc) ~ # No LFC in output (DEXSeq)
                  dplyr::summarise(., pvalue = aggregation::lancaster(pvalue,
                                                                      .data[[weights]])),
                !is.null(lfc) ~ # With LFC in output (DESeq2)
                  dplyr::summarise(., pvalue = aggregation::lancaster(pvalue,
                                                                      .data[[weights]]),
                                   lfc = weighted.mean(.data[[lfc]], .data[[weights]]))) %>% 
    # Remove zeros again to prevent downstream issues
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) 
  return(res)
}
