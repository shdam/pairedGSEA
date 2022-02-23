
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
  
  ### Check that results exists
  deseq2file <- paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS")
  if(!file.exists(deseq2file)) stop(paste0("File:", deseq2file, " does not exists.\\n","Please run experiment before analysing. See ?runExperiment"))
  
  
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  
  
  message("Analysing ", row$study, " ", experimentTitle)
  # if(FALSE){
  # Load results
  message("Loading results")
  
  res <- readRDS(paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
  dxr <- readRDS(paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
  
  # p value aggregation
  
  message("Aggregating p values")
  
  dxr_agg <- perGenePValue(dxr, gene = "groupID", weights = "exonBaseMean", lfc = "log2FC_baseline_vs_condition", type = "dexseq")# %>% 
    # filter(pvalue < 0.1)
  res_agg <- perGenePValue(res, gene = "gene", weights = "baseMean", lfc = "log2FC", type = "deseq2")# %>% 
    # filter(pvalue < 0.1)
  
  (comb <- dplyr::full_join(res_agg, dxr_agg, by = "ensembl_gene", suffix = c("_deseq2", "_dexseq")))
  
  
  message("Storing aggregation result")
  saveRDS(comb, paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
  # }
  # comb <- readRDS(paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
  if(FALSE){
  message("Gene set enrichment analysis")
  
  ### Defining stats
  resStats <- comb %>% 
    dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_deseq2) * sign(lfc_deseq2)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene)

  
  dxrStats <- comb %>% 
    dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_dexseq)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene) 
  dxrStats2 <- comb %>% 
    dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_dexseq) * sign(lfc_dexseq)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene) 
  
  
  resStatss <- as.vector( scale( rev(seq_along( resStats)), c=T,s=T))
  names(resStatss) <- names(resStats)
  
  dxrStatss <- as.vector( scale( rank( dxrStats), c=T,s=T))
  names(dxrStatss) <- names(dxrStats)
  dxrStatss2 <- as.vector( scale( rank( dxrStats2), c=T,s=T))
  names(dxrStatss2) <- names(dxrStats2)
  
  ### Run fgsea
  message("Running fgsea on DESeq2 results")
  fgseaRes_oristd <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = resStats,
                                     scoreType = "std",
                                     eps = 10e-320
  )

  # saveRDS(fgseaRes, paste0("results/", dataname, "_fgseaRes_", experimentTitle, ".RDS"))
  # fgseaRes_sig <- fgseaRes %>% filter(padj<0.05)
  message("Running fgsea on DEXSeq results")
  fgseaDxr_pos <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = dxrStatss,
                                     scoreType = "pos",
                                     eps = 10e-320
  )
  fgseaDxr_std <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = dxrStatss2,
                                     scoreType = "std",
                                     eps = 10e-320
  )
  # saveRDS(fgseaDxr, paste0("results/", dataname, "_fgseaDxr_", experimentTitle, ".RDS"))
  # fgseaDxr_sig <- fgseaDxr %>% filter(padj<0.05)
  fgseaDxr22 <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                      stats = dxrStatss,
                                      scoreType = "pos",
                                      eps = 10e-320
  )
  # fgseaDxr_sig2 <- fgseaDxr2 %>% filter(padj<0.05)
  saveRDS(fgseaDxr2, paste0("results/", dataname, "_fgseaDxr2_", experimentTitle, ".RDS"))
  message("fgsea results are stored in the results folder. Look for '*_fgseaRes_*' and '*_fgseaDxr_*'")} # fgsea
  
  if(TRUE){ # fora
    message("Running fora")
    ### Significant genes
    resGenes <- comb %>% 
      filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>%
      mutate(padj = p.adjust(pvalue_deseq2, "fdr")) %>% 
      filter(padj < 0.05) %>% 
      arrange(padj)
    dxrGenes <- comb %>% 
      filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>%
      mutate(padj = p.adjust(pvalue_dexseq, "fdr")) %>% 
      filter(padj < 0.05) %>% 
      arrange(padj)
    
    forares <- fgsea::fora(gene_sets, genes = resGenes$ensembl_gene, 
                           universe = unique(comb$ensembl_gene), minSize = 25)
    foradxr <- fgsea::fora(gene_sets, genes = dxrGenes$ensembl_gene,
                           universe = unique(comb$ensembl_gene), minSize = 25)
    foraresdxr <- fgsea::fora(gene_sets, genes = intersect(dxrGenes$ensembl_gene, resGenes$ensembl_gene,),
                              universe = unique(comb$ensembl_gene), minSize = 25)
    message("Storing fora results")
    saveRDS(forares, paste0("results/", dataname, "_forares_", experimentTitle, ".RDS"))
    saveRDS(foradxr, paste0("results/", dataname, "_foradxr_", experimentTitle, ".RDS"))
    saveRDS(foraresdxr, paste0("results/", dataname, "_foraresdxr_", experimentTitle, ".RDS"))
  }
  
  
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
                           lfc = "log2FC",
                           type = "deseq2"){
  stopifnot(typeof(df) == "list")
  stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))
  type <- stringr::str_to_lower(type)
  res <- df %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::rename(pvalue = .data[[p]],
                  ensembl_gene = .data[[gene]],
                  lfc = .data[[lfc]],
                  weights = .data[[weights]]) %>% 
    # Prevent warning from Lancaster
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) %>% 
    dplyr::group_by(ensembl_gene) %>% 
    # p value aggregation
    purrr::when(type == "dexseq" ~ 
                  dplyr::summarise(.,
                                   lfc = lfc[which(pvalue == min(pvalue))],
                                   # stat = stat[which(pvalue == min(pvalue))],
                                   pvalue = aggregation::lancaster(pvalue, weights)),
                type == "deseq2" ~ 
                  dplyr::summarise(.,
                                   lfc = weighted.mean(lfc, weights),
                                   # stat2 = weighted.mean(stat, weights),
                                   # stat = stat[which(pvalue == min(pvalue))],
                                   pvalue = aggregation::lancaster(pvalue, weights))) %>% 
    # Remove zeros again to prevent downstream issues
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) 
  return(res)
  }
