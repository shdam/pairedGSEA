
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

    message("Adding ", row$study, " ", experimentTitle)
    
    ### Load results
    fgsea_deseq <- readRDS(paste0("results/", experiment_title, "_fgsea_deseq.RDS"))
    fgsea_dexseq <- readRDS(paste0("results/", experiment_title, "_fgsea_dexseq.RDS"))
    fgsea_paired <- readRDS(paste0("results/", experiment_title, "_fgsea_paired.RDS"))
    
    # fgseatot <- fgseatot <-
    #   fgseaRes %>% 
    #   
    
    ### Significant gene sets
    pathRes <- fgseaRes %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathDxr <- fgseaDxr2 %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathDxr2 <- fgseaDxr3 %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    
    overlap <- union(pathRes, pathDxr) %>% union(pathDxr2) %>% 
      mutate(deseq2 = pathway %in% pathRes$pathway,
             dexseq = pathway %in% pathDxr$pathway,
             dexseq2 = pathway %in% pathDxr2$pathway,
             overlap = deseq2 + dexseq + dexseq2,
             experiment = paste(row$study, experimentTitle)
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
concatenate_fora <- function(experiments){
  
  fora_all <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove(".csv") 
    ### Define experiment details
    comparison <- row$`comparison (baseline_v_condition)` %>% pairedGSEA:::check_comparison()
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    
    message("Adding ", experiment_title)
    
    ### Load results
    fora_deseq <- readRDS(paste0("results/", experiment_title, "_fora_deseq.RDS"))
    fora_dexseq <- readRDS(paste0("results/", experiment_title, "_fora_dexseq.RDS"))
    fora_paired <- readRDS(paste0("results/", experiment_title, "_fora_paired.RDS"))
    
    ### Join results
    fora_joined <- fora_deseq %>% 
      dplyr::full_join(fora_dexseq, by = c("pathway", "size_geneset", "size_universe"), suffix = c("_deseq", "_dexseq")) %>% 
      dplyr::full_join(fora_paired, by = c("pathway", "size_geneset", "size_universe")) %>% 
      dplyr::rename(padj_paired = padj,
                    pval_paired = pval,
                    overlap_paired = overlap,
                    size_genes_paired = size_genes,
                    odds_ratio_paired = odds_ratio,
                    enrichment_score_paired = enrichment_score) %>% 
      dplyr::select(-dplyr::starts_with("overlapGenes")) %>% 
      dplyr::mutate(experiment = experiment_title)
    
    fora_all <- fora_all %>% 
      dplyr::bind_rows(fora_joined)
  }
  saveRDS(fora_all, "results/fora_all.RDS")
  return(fora_all)
}
# deseq2 + dexseq ----
concatRes <- function(experiments){
  
  concatResults <- tibble::tibble()
  
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
    message("Adding ", row$study, " ", experimentTitle)
    res <- readRDS(paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
    dxr <- readRDS(paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
    
    comb <- res %>% 
      dplyr::left_join(dxr, by = c("transcript" = "featureID"), suffix = c("_deseq2", "_dexseq")) %>% 
      dplyr::mutate(experiment = paste(row$study, experimentTitle)) %>% 
      dplyr::select(-dplyr::starts_with("count"))
    
    
    concatResults <- concatResults %>% 
      dplyr::bind_rows(comb)
  }
  saveRDS(concatResults, "results/concatResults.RDS")
  return(concatResults)
}


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
    comparison <- row$`comparison (baseline_v_condition)` %>% pairedGSEA:::check_comparison()
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    ### Check that results exists
    message("Adding ", experiment_title)
    comb <- readRDS(paste0("results/", experiment_title, "_aggregated_pvals.RDS"))
    comb$experiment <- experiment_title
    
    
    concatenated_genes <- concatenated_genes %>% 
      dplyr::bind_rows(comb)
  }
  concatenated_genes <- concatenated_genes %>% 
    dplyr::mutate(padj_deseq = p.adjust(pvalue_deseq, "fdr"),
                  padj_dexseq = p.adjust(pvalue_dexseq, "fdr"))
  saveRDS(concatenated_genes, "results/concatenated_genes.RDS")
  return(concatenated_genes)
}


# Concat SVA level ----
concatSVA <- function(experiments){
  
  concatsva <- tibble::tibble()
  
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
    message("Adding ", row$study, " ", experimentTitle)
    dds <- readRDS(paste0("results/", dataname, "_dds_", experimentTitle, ".RDS"))
    des <- tibble::tribble(~experiment, ~design,
                           paste(dataname, experimentTitle), DESeq2::design(dds))
    
    
    concatsva <- concatsva %>% 
      dplyr::bind_rows(des)
  }
  saveRDS(concatsva, "results/concatSVA.RDS")
  return(concatsva)
}
