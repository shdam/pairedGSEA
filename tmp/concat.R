
# fgsea ----
concatFgseaResults <- function(experiments){
  
  concatFgsea <- tibble::tibble()
  fgseatot <- tibble::tibble()
  
  for(row in 1:nrow(experiments)){
    row <- experiments[row, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    comparison <- row$`comparison (baseline_v_condition)` %>% pairedGSEA:::check_comparison()
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)

    message("Adding ", row$study, " ", experiment_title)
    
    ### Load results
    fgsea_deseq <- readRDS(paste0("results/", experiment_title, "_fgsea_deseq.RDS"))
    fgsea_dexseq <- readRDS(paste0("results/", experiment_title, "_fgsea_dexseq.RDS"))
    fgsea_dexseqlfc <- readRDS(paste0("results/", experiment_title, "_fgsea_dexseqlfc.RDS"))
    
    # fgseatot <- fgseatot <-
    #   fgseaRes %>% 
    #   
    
    ### Significant gene sets
    pathRes <- fgsea_deseq %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathDxr <- fgsea_dexseq %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathDxr2 <- fgsea_dexseqlfc %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    
    overlap <- union(pathRes, pathDxr) %>% union(pathDxr2) %>% 
      mutate(deseq = pathway %in% pathRes$pathway,
             dexseq = pathway %in% pathDxr$pathway,
             dexseqlfc = pathway %in% pathDxr2$pathway,
             overlap = deseq + dexseq + dexseqlfc,
             experiment = experiment_title
      )
    concatFgsea <- concatFgsea %>% 
      dplyr::bind_rows(overlap)
  }
  saveRDS(concatFgsea, "results/concatFgsea.RDS")
  return(concatFgsea)
}

# fora ----
concatForaResults <- function(experiments){
  
  concatFora <- tibble::tibble()
  foratot <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
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
    
    message("Adding ", row$study, " ", experimentTitle)
    
    ### Load results
    forares <- readRDS(paste0("results/", dataname, "_forares_", experimentTitle, ".RDS"))
    foradxr <- readRDS(paste0("results/", dataname, "_foradxr_", experimentTitle, ".RDS"))
    foraresdxr <- readRDS(paste0("results/", dataname, "_foraresdxr_", experimentTitle, ".RDS"))
    
    comb <- forares %>% 
      dplyr::rename(padj_deseq2 = padj,
             pval_deseq2 = pval,
             overlap_deseq2 = overlap) %>% 
      dplyr::full_join(foradxr, by = "pathway") %>% 
      dplyr::rename(padj_dexseq = padj,
             pval_dexseq = pval,
             overlap_dexseq = overlap) %>% 
      dplyr::full_join(foraresdxr, by = "pathway") %>% 
      dplyr::rename(padj_decombined = padj,
             pval_decombined = pval,
             overlap_decombined = overlap) %>% 
      dplyr::select(pathway, starts_with("padj_"), starts_with("pval_"), starts_with("overlap_"), size) %>% 
      dplyr::mutate(experiment = paste(row$study, experimentTitle))
    foratot <- foratot %>% 
      dplyr::bind_rows(comb)
    
    
    ### Significant gene sets
    pathres <- forares %>% dplyr::filter(padj < 0.05) %>% dplyr::select(pathway) %>% tibble::as_tibble()
    pathdxr <- foradxr %>% dplyr::filter(padj < 0.05) %>% dplyr::select(pathway) %>% tibble::as_tibble()
    pathresdxr <- foraresdxr %>% dplyr::filter(padj < 0.05) %>% dplyr::select(pathway) %>% tibble::as_tibble()
    
    overlap <- dplyr::union(pathres, pathdxr) %>% dplyr::union(pathresdxr) %>% 
      dplyr::mutate(deseq2 = pathway %in% pathres$pathway,
             dexseq = pathway %in% pathdxr$pathway,
             decombined = pathway %in% pathresdxr$pathway,
             overlap = deseq2 + dexseq + decombined,
             experiment = paste(row$study, experimentTitle)
      )
    concatFora <- concatFora %>% 
      dplyr::bind_rows(overlap)
  }
  saveRDS(foratot, "results/foratot.RDS")
  saveRDS(concatFora, "results/concatFora.RDS")
  return(concatFora)
}
concatenate_ora <- function(experiments){
  
  ora_all <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove(".csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    
    message("Adding ", experiment_title)
    
    ### Load results
    # fora_deseq <- readRDS(paste0("results/", experiment_title, "_fora_deseq.RDS"))
    # fora_dexseq <- readRDS(paste0("results/", experiment_title, "_fora_dexseq.RDS"))
    ora_joined <- readRDS(paste0("results/", experiment_title, "_ora.RDS"))
    
    ### Join results
    # fora_joined <- fora_deseq %>% 
    #   dplyr::full_join(fora_dexseq, by = c("pathway", "size_geneset", "size_universe"), suffix = c("_deseq", "_dexseq")) %>% 
    #   dplyr::full_join(fora_paired, by = c("pathway", "size_geneset", "size_universe")) %>% 
    #   dplyr::rename(padj_paired = padj,
    #                 pval_paired = pval,
    #                 overlap_paired = overlap,
    #                 size_genes_paired = size_genes,
    #                 odds_ratio_paired = odds_ratio,
    #                 enrichment_score_paired = enrichment_score) %>% 
    #   dplyr::select(-dplyr::starts_with("overlapGenes")) %>% 
    #   dplyr::mutate(experiment = experiment_title)
    
    ora_all <- ora_all %>% 
      dplyr::bind_rows(ora_joined)
  }
  saveRDS(ora_all, "results/ora_all.RDS")
  return(ora_all)
}
# concatenated_ora <- concatenate_ora(experiments)

# deseq2 + dexseq ----
concatenate_results <- function(experiments){
  
  concatenated_results <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    ### Check that results exists
    message("Adding ", row$study, " ", experiment_title)
    results_deseq <- readRDS(paste0("results/", experiment_title, "_deseqres.RDS"))
    results_dexseq <- readRDS(paste0("results/", experiment_title, "_dexseqres.RDS"))
    
    results_dexseq <- results_dexseq %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(lfc_dexseq = log2fold_C_B) %>% 
    results_deseq <- results_deseq %>%
      tibble::as_tibble(rownames = "gene_tx") %>% 
      dplyr::rename(lfc_deseq = log2FoldChange) %>% 
      tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") 
    
    genes_to_keep_dexseq <- results_dexseq %>%
      dplyr::group_by(groupID) %>% 
      dplyr::summarise(n_significant = sum(padj < 0.05, na.rm = TRUE)) %>% 
      dplyr::filter(n_significant > 1) %>% 
      dplyr::pull(groupID)
    genes_to_keep_deseq <- results_deseq %>%
      # tibble::as_tibble(rownames = "gene_tx") %>% 
      # tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
      # dplyr::rename(log2FC_deseq = log2FoldChange) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::summarise(n_significant = sum(padj < 0.05, na.rm = TRUE)) %>% 
      dplyr::filter(n_significant > 1) %>% 
      dplyr::pull(gene)
    
    genes_to_keep <- union(genes_to_keep_dexseq, genes_to_keep_deseq)
    
    results_joined <- results_deseq %>% 
      dplyr::left_join(results_dexseq, by = c("transcript" = "featureID"), suffix = c("_deseq", "_dexseq")) %>% 
      dplyr::filter(gene %in% genes_to_keep) %>% 
      dplyr::mutate(experiment = experiment_title) %>% 
      dplyr::select(-dplyr::starts_with("count"))
    
    
    concatenated_results <- concatenated_results %>% 
      dplyr::bind_rows(results_joined)
  }
  saveRDS(concatenated_results, "results/concatenated_results.RDS")
  return(concatenated_results)
}
# concatenated_results <- concatenate_results(experiments)

# Concat gene level ----
concatenate_genes <- function(experiments){
  
  concatenated_genes <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    ### Loading results
    message("Adding ", experiment_title)
    aggregated_pvals <- readRDS(paste0("results/", experiment_title, "_aggregated_pvals.RDS"))
    aggregated_pvals <- aggregated_pvals %>% 
      dplyr::mutate(padj_deseq = p.adjust(pvalue_deseq, "fdr"),
                    padj_dexseq = p.adjust(pvalue_dexseq, "fdr"),
                    experiment = experiment_title)
    
    
    concatenated_genes <- concatenated_genes %>% 
      dplyr::bind_rows(aggregated_pvals)
  }

  saveRDS(concatenated_genes, "results/concatenated_genes.RDS")
  return(concatenated_genes)
}

concatenate_no_sva_genes <- function(experiments){
  
  concatenated_no_sva_genes <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`, "_no_sva")
    ### Loading results
    message("Adding ", experiment_title)
    aggregated_pvals <- readRDS(paste0("results/", experiment_title, "_aggregated_pvals.RDS"))
    aggregated_pvals <- aggregated_pvals %>% 
      dplyr::mutate(experiment = stringr::str_remove(experiment_title, "_no_sva"))
    
    
    concatenated_no_sva_genes <- concatenated_no_sva_genes %>% 
      dplyr::bind_rows(aggregated_pvals)
  }
  
  saveRDS(concatenated_no_sva_genes, "results/concatenated_no_sva_genes.RDS")
  return(concatenated_no_sva_genes)
}
# concatenated_no_sva_genes <- concatenate_no_sva_genes(experiments)

# Concat SVA level ----
concatenate_design <- function(experiments){
  
  concatenated_design <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    ### Check that results exists
    message("Adding ", experiment_title)
    dds <- readRDS(paste0("results/", experiment_title, "_dds.RDS"))
    des <- tibble::tribble(~experiment, ~design,
                           experiment_title, DESeq2::design(dds))
    
    
    concatenated_design <- concatenated_design %>% 
      dplyr::bind_rows(des)
  }
  saveRDS(concatenated_design, "results/concatenated_design.RDS")
  return(concatenated_design)
}
# concatenated_design <- concatenate_design(experiments)

# deseq_genes_no_sva <- deseq_no_sva(experiments, archs4db, gtf)


concatenate_sva_correlation <- function(experiments){
  
  concatenated_correlation <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    ### Check that results exists
    message("Adding ", experiment_title)
    dds <- readRDS(paste0("results/", experiment_title, "_dds.RDS"))
    md <- colData(dds) %>% dplyr::as_tibble()
    correlation <- apply(select(md, starts_with("sv")), 2, function(x) cor(as.numeric(md$group_nr), x))
    correlation$experiment <- experiment_title
    
    concatenated_correlation <- concatenated_correlation %>% 
      dplyr::bind_rows(correlation)
  }
  saveRDS(concatenated_correlation, "results/concatenated_sva_correlation.RDS")
  return(concatenated_correlation)
}
# cors <- concatenaded_correlation %>% pivot_longer(cols = starts_with("sv"), names_to = "sv", values_to = "correlation")



# Store results for paper ----
store_experiment_results <- function(experiments){
  check_make_dir("results/paper")
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    
    message("Storing ", experiment_title)
    metadata <- readxl::read_xlsx(md_file)
    
    aggregated_pvals <- readRDS(paste0("results/", experiment_title, "_aggregated_pvals.RDS"))
    ora <- readRDS(paste0("results/", experiment_title, "_ora.RDS"))
    
    store_result(list("metadata" = metadata,
                      "genes" = aggregated_pvals,
                      "gene_sets" = ora),
                 file = paste0("paper/", experiment_title, ".RDS"),
                 analysis = experiment_title)
  }
  message("Done!\n")
}
# store_experiment_results(experiments)