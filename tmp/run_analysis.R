




#' Analyse experiments
#' @noRd
run_analysis <- function(row){
  
  
  if(typeof(row) == "character"){ # Convert apply-made row to tibble
    row <- tibble::as_tibble(row, rownames = "names") %>% 
      tidyr::pivot_wider(values_from = value, names_from = names)
  }
  
  
  ### Load metadata
  md_file <- row$filename
  data_name <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove(".csv") 
  
  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)` %>% pairedGSEA:::check_comparison()
  experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
  
  ### Check that results exists
  # deseq2file <- paste0("results/", experimentTitle, "_deseq2res.RDS")
  # if(!file.exists(deseq2file)) stop(paste0("File:", experimentTitle, " does not exists.\\n","Please run experiment before analysing. See ?runExperiment"))
  
  
  message("Analysing ", experiment_title)
  if(TRUE){
    # Load results
    message("Loading results")
    
    res <- readRDS(paste0("results/", experiment_title, "_deseq2res.RDS"))
    dxr <- readRDS(paste0("results/", experiment_title, "_dexseqres.RDS"))
    
    # p value aggregation
    
    message("Aggregating p values")
    
    dxr_agg <- per_gene_pvalue(dxr, gene = "groupID", weights = "exonBaseMean", lfc = "log2FC_dexseq", type = "dexseq")
    res_agg <- per_gene_pvalue(res, gene = "gene", weights = "baseMean", lfc = "log2FC_deseq", type = "deseq")
    
    (aggregated_pvals <- dplyr::full_join(res_agg, dxr_agg, by = "ensembl_gene", suffix = c("_deseq", "_dexseq")))
    
    
    # message("Storing aggregation result")
    pairedGSEA:::store_result(aggregated_pvals, paste0(experiment_title, "_aggregated_pvals.RDS"), "gene aggregation")
  } else{
    comb <- readRDS(paste0(experiment_title, "_aggregated_pvals.RDS"))
  }
  
  if(FALSE){
    message("Gene set enrichment analysis")
    
    ### Defining stats
    stats_deseq <- aggregated_pvals %>% 
      dplyr::filter(!is.na(pvalue_deseq) & !is.na(ensembl_gene)) %>% 
      dplyr::mutate(pvalue = -log10(pvalue_deseq) * sign(lfc_deseq)) %>% 
      dplyr::pull(pvalue, name = ensembl_gene)
    
    
    stats_dexseq <- aggregated_pvals %>% 
      dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
      dplyr::mutate(pvalue = -log10(pvalue_dexseq)) %>% 
      dplyr::pull(pvalue, name = ensembl_gene) 
    stats2_dexseq <- aggregated_pvals %>% 
      dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
      dplyr::mutate(pvalue = -log10(pvalue_dexseq) * sign(lfc_dexseq)) %>% 
      dplyr::pull(pvalue, name = ensembl_gene) 
    
    
    ### Run fgsea
    message("Running fgsea on DESeq2 results")
    fgsea_deseq <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                          stats = stats_deseq,
                                          nproc = 10,
                                          scoreType = "std",
                                          eps = 10e-320,
                                          minSize = 25
    )
    pairedGSEA:::store_result(fgsea_deseq, paste0(experiment_title, "_fgsea_deseq.RDS"), "fgsea on deseq2 results")
    
    message("Running fgsea on DEXSeq results")
    fgsea_dexseq <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                           stats = stats_dexseq,
                                           nproc = 10,
                                           scoreType = "std",
                                           eps = 10e-320,
                                           minSize = 25
    )
    fgsea_lfc_dexseq <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                               stats = stats2_dexseq,
                                               nproc = 10,
                                               scoreType = "std",
                                               eps = 10e-320,
                                               minSize = 25
    )
    pairedGSEA:::store_result(fgsea_dexseq, paste0(experiment_title, "_fgsea_dexseq.RDS"), "fgsea on dexseq results")
    pairedGSEA:::store_result(fgsea_lfc_dexseq, paste0(experiment_title, "_fgsea_dexseqlfc.RDS"), "fgsea on dexseq results using LFC")
    
    message("fgsea results are stored in the results folder. Look for '*_fgsea_deseq.RDS' and '*_fgsea_dexseq*.RDS'")
  } # fgsea
  
  if(TRUE){ # fora
    message("Running fora")
    ### Significant genes
    genes_deseq <- aggregated_pvals %>% 
      dplyr::filter(!is.na(pvalue_deseq) & !is.na(ensembl_gene)) %>%
      dplyr::mutate(padj = stats::p.adjust(pvalue_deseq, "fdr")) %>% 
      dplyr::filter(padj < 0.05) %>% 
      dplyr::arrange(padj)
    genes_dexseq <- aggregated_pvals %>% 
      dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>%
      dplyr::mutate(padj = stats::p.adjust(pvalue_dexseq, "fdr")) %>% 
      dplyr::filter(padj < 0.05) %>% 
      dplyr::arrange(padj)
    genes_paired <- genes_deseq %>% 
      dplyr::full_join(genes_dexseq, by = "ensembl_gene")
    
    fora_deseq <- fgsea::fora(gene_sets, genes = genes_deseq$ensembl_gene, 
                              universe = unique(aggregated_pvals$ensembl_gene), minSize = 25) %>% 
      dplyr::rename(size_geneset = size) %>% 
      dplyr::mutate(size_genes = length(genes_deseq$ensembl_gene),
                    size_universe = length(unique(aggregated_pvals$ensembl_gene)),
                    odds_ratio = (overlap / size_genes) / (size_geneset / size_universe),
                    enrichment_score = log2(odds_ratio)) 
    
    fora_dexseq <- fgsea::fora(gene_sets, genes = genes_dexseq$ensembl_gene,
                               universe = unique(aggregated_pvals$ensembl_gene), minSize = 25) %>% 
      dplyr::rename(size_geneset = size) %>% 
      dplyr::mutate(size_genes = length(genes_dexseq$ensembl_gene),
                    size_universe = length(unique(aggregated_pvals$ensembl_gene)),
                    odds_ratio = (overlap / size_genes) / (size_geneset / size_universe),
                    enrichment_score = log2(odds_ratio))
    
    fora_paired <- fgsea::fora(gene_sets, genes = genes_paired$ensembl_gene,
                               universe = unique(aggregated_pvals$ensembl_gene), minSize = 25) %>% 
      dplyr::rename(size_geneset = size) %>% 
      dplyr::mutate(size_genes = length(genes_paired$ensembl_gene),
                    size_universe = length(unique(aggregated_pvals$ensembl_gene)),
                    odds_ratio = (overlap / size_genes) / (size_geneset / size_universe),
                    enrichment_score = log2(odds_ratio))
    
    message("Storing fora results")
    pairedGSEA:::store_result(fora_deseq, paste0(experiment_title, "_fora_deseq.RDS"), "fora on deseq2 results")
    pairedGSEA:::store_result(fora_dexseq, paste0(experiment_title, "_fora_dexseq.RDS"), "fora on dexseq results")
    pairedGSEA:::store_result(fora_paired, paste0(experiment_title, "_fora_paired.RDS"), "fora on paired results")
  }
  
  
}

#' Per gene p-value aggregation
#' @noRd
perGenePValue <- function (df,
                           gene = "gene",
                           p = "pvalue",
                           weights = "baseMean",
                           lfc = "log2FC",
                           type = "deseq2"){
  stopifnot("data.frame" %in% class(df))
  stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))
  type <- tolower(type)
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
                                   lfc = lfc[which(pvalue == min(pvalue))][[1]],
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
#' Per gene p-value aggregation
#' @noRd
per_gene_pvalue <- function (df,
                             gene = "gene",
                             p = "pvalue",
                             weights = "baseMean",
                             lfc = "log2FC",
                             type = "deseq"){
  stopifnot("data.frame" %in% class(df))
  stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))
  type <- tolower(type)
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
                                   lfc = lfc[which(pvalue == min(pvalue))][[1]],
                                   # stat = stat[which(pvalue == min(pvalue))],
                                   pvalue = aggregation::lancaster(pvalue, weights)),
                type == "deseq" ~ 
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



#' Run SVA and export dds
#' @noRd
getDDS <- function(row, archs4db = NULL, tx_count = NULL, group_col = "group_nr", tpm = TRUE, prefilter = 10, parallel = TRUE){
  
  if(typeof(row) == "character"){ # Convert apply-made row to tibble
    row <- tibble::as_tibble(row, rownames = "names") %>% 
      tidyr::pivot_wider(values_from = value, names_from = names)
  }
  
  
  message("Running on ", row$study)
  
  ### Load metadata
  md_file <- row$filename
  dataname <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove(".csv") 
  
  ### Define tpm file
  if(tpm) tpm <- stringr::str_replace(archs4db, "counts", "tpm")
  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)` %>% pairedGSEA:::check_comparison()
  experiment_title <- row$`comparison_title (empty_if_not_okay)`
  
  if(!is.null(archs4db)) tx_count <- prepare_tx_count(md_file,
                                                      group_col,
                                                      comparison,
                                                      archs4db)
  # Check parallel
  if(parallel) pairedGSEA:::check_missing_package("BiocParallel", "Bioc")
  
  # Load metadata
  message("Preparing metadata")
  metadata <- pairedGSEA:::prepare_metadata(metadata, group_col, comparison)
  
  # Subsample in case metadata 
  metadata <- metadata[metadata[[sample_col]] %in% colnames(tx_count), ]
  # Ensure rows in metadata matches columns in the count matrix
  tx_count <- tx_count[, metadata[[sample_col]]]
  
  ### Prefiltering
  if(prefilter) tx_count <- pairedGSEA:::pre_filter(tx_count, prefilter)
  
  if(!quiet) message("Converting count matrix to DESeqDataSet")
  # Create DDS from count matrix
  dds <- pairedGSEA:::convert_matrix_to_dds(tx_count, metadata, group_col)
  
  dds <- pairedGSEA:::run_sva(dds, group_col, quiet = quiet)
  
  
  saveRDS(dds, paste0("results/", dataname, "_dds_", experimentTitle, ".RDS"))
  
  return(dds)
}



#' Analyse experiments
#' @noRd
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
  comparison <- row$`comparison (baseline_v_condition)` #%>% pairedGSEA:::check_comparison()
  experimentTitle <- row$`comparison_title (empty_if_not_okay)`
  
  ### Check that results exists
  # deseq2file <- paste0("results/", experimentTitle, "_deseq2res.RDS")
  # if(!file.exists(deseq2file)) stop(paste0("File:", experimentTitle, " does not exists.\\n","Please run experiment before analysing. See ?runExperiment"))
  
  
  message("Analysing ", row$study, " ", experimentTitle)
  if(FALSE){
    # Load results
    message("Loading results")
    
    res <- readRDS(paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
    dxr <- readRDS(paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
    
    # p value aggregation
    
    message("Aggregating p values")
    
    dxr_agg <- perGenePValue(dxr, gene = "groupID", weights = "exonBaseMean", lfc = "log2FC_baseline_vs_condition", type = "dexseq")
    res_agg <- perGenePValue(res, gene = "gene", weights = "baseMean", lfc = "log2FC", type = "deseq2")
    
    (comb <- dplyr::full_join(res_agg, dxr_agg, by = "ensembl_gene", suffix = c("_deseq2", "_dexseq")))
    
    
    message("Storing aggregation result")
    saveRDS(comb, paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
  } else{
    comb <- readRDS(paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
  }
  
  if(TRUE){
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
    
    # 
    # resStatss <- as.vector( scale( rev(seq_along( resStats)), c=T,s=T))
    # names(resStatss) <- names(resStats)
    # 
    # dxrStatss <- as.vector( scale( rank( dxrStats), c=T,s=T))
    # names(dxrStatss) <- names(dxrStats)
    # dxrStatss2 <- as.vector( scale( rank( dxrStats2), c=T,s=T))
    # names(dxrStatss2) <- names(dxrStats2)
    
    if(FALSE){
      fgseaDxr3_std <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                              stats = dxrStats,
                                              nproc = 10,
                                              scoreType = "std",
                                              eps = 10e-320,
                                              minSize = 25
      )
      saveRDS(fgseaDxr3_std, paste0("results/", dataname, "_fgseaDxr3_", experimentTitle, ".RDS"))
    } else{
      
      
      
      ### Run fgsea
      message("Running fgsea on DESeq2 results")
      fgseaRes_oristd <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                                stats = resStats,
                                                nproc = 10,
                                                scoreType = "std",
                                                eps = 10e-320,
                                                minSize = 25
      )
      
      saveRDS(fgseaRes_oristd, paste0("results/", dataname, "_fgseaRes_", experimentTitle, ".RDS"))
      # fgseaRes_sig <- fgseaRes %>% filter(padj<0.05)
      message("Running fgsea on DEXSeq results")
      fgseaDxr_pos <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                             stats = dxrStats,
                                             nproc = 10,
                                             scoreType = "pos",
                                             eps = 10e-320,
                                             minSize = 25
      )
      fgseaDxr_std <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                             stats = dxrStats2,
                                             nproc = 10,
                                             scoreType = "std",
                                             eps = 10e-320,
                                             minSize = 25
      )
      fgseaDxr3_std <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                              stats = dxrStats,
                                              nproc = 10,
                                              scoreType = "std",
                                              eps = 10e-320,
                                              minSize = 25
      )
      saveRDS(fgseaDxr_pos, paste0("results/", dataname, "_fgseaDxr_", experimentTitle, ".RDS"))
      saveRDS(fgseaDxr_std, paste0("results/", dataname, "_fgseaDxr2_", experimentTitle, ".RDS"))
      saveRDS(fgseaDxr3_std, paste0("results/", dataname, "_fgseaDxr3_", experimentTitle, ".RDS"))
    }
    # fgseaDxr_sig <- fgseaDxr %>% filter(padj<0.05)
    # fgseaDxr22 <- fgsea::fgseaMultilevel(pathways = gene_sets,
    #                                     stats = dxrStatss,
    #                                     scoreType = "pos",
    #                                     eps = 10e-320
    # )
    # fgseaDxr_sig2 <- fgseaDxr2 %>% filter(padj<0.05)
    # saveRDS(fgseaDxr2, paste0("results/", dataname, "_fgseaDxr2_", experimentTitle, ".RDS"))
    message("fgsea results are stored in the results folder. Look for '*_fgseaRes_*' and '*_fgseaDxr*'")
  } # fgsea
  
  if(FALSE){ # fora
    message("Running fora")
    ### Significant genes
    resGenes <- comb %>% 
      filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>%
      mutate(padj = stats::p.adjust(pvalue_deseq2, "fdr")) %>% 
      filter(padj < 0.05) %>% 
      arrange(padj)
    dxrGenes <- comb %>% 
      filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>%
      mutate(padj = stats::p.adjust(pvalue_dexseq, "fdr")) %>% 
      filter(padj < 0.05) %>% 
      arrange(padj)
    resdxrGenes <- resGenes %>% 
      full_join(dxrGenes, by = "ensembl_gene")
    
    forares <- fgsea::fora(gene_sets, genes = resGenes$ensembl_gene, 
                           universe = unique(comb$ensembl_gene), minSize = 25)
    foradxr <- fgsea::fora(gene_sets, genes = dxrGenes$ensembl_gene,
                           universe = unique(comb$ensembl_gene), minSize = 25)
    foraresdxr <- fgsea::fora(gene_sets, genes = resdxrGenes$ensembl_gene,
                              universe = unique(comb$ensembl_gene), minSize = 25)
    message("Storing fora results")
    saveRDS(forares, paste0("results/", dataname, "_forares_", experimentTitle, ".RDS"))
    saveRDS(foradxr, paste0("results/", dataname, "_foradxr_", experimentTitle, ".RDS"))
    saveRDS(foraresdxr, paste0("results/", dataname, "_foraresdxr_", experimentTitle, ".RDS"))
  }
  
  
}
