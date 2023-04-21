## Concat all ----

if(FALSE){
  
  experiments <- combine_experiments()
  
  concatenate_ora(experiments)
  concatenate_results(experiments)
  concatenate_genes(experiments)
  concatenate_no_sva_genes(experiments)
  concatenate_design(experiments)
  concatenate_sva_correlation(experiments)
  
  
  experiments_limma <- combine_experiments(limma = TRUE)
  
  concatenate_ora(experiments_limma, suffix = "_limma_all")
  concatenate_results(experiments_limma, suffix = "_limma_results", limma = TRUE)
  concatenate_genes(experiments_limma, suffix = "_limma_genes", limma = TRUE)
  concatenate_no_sva_genes(experiments_limma, suffix = "_limma_no_sva_genes")
  # concatenate_design(experiments_limma)
  # concatenate_sva_correlation(experiments_limma)
  
}


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

concatenate_ora <- function(experiments, suffix = "_all"){
  
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
  saveRDS(ora_all, paste0("results/ora", suffix,".RDS"))
  return(ora_all)
}
# concatenated_ora <- concatenate_ora(experiments)

# deseq2 + dexseq ----
concatenate_results <- function(experiments, suffix = "_results", new = FALSE, limma = FALSE){
  
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
    if(new | limma){
      suffix_expression <- "_expression"
      suffix_splicing <- "_splicing"
    } else{
      suffix_expression <- "_deseqres"
      suffix_splicing <- "_dexseqres"
    }
    ### Check that results exists
    message("Adding ", row$study, " ", experiment_title)
    if(limma){
      results <- readRDS(paste0("results/", experiment_title, suffix, ".RDS"))
      results_expression <- results$expression
      results_splicing <- results$splicing
      rm(results)
    } else{
      results_expression <- readRDS(paste0("results/", experiment_title, suffix_expression, ".RDS"))
      results_splicing <- readRDS(paste0("results/", experiment_title, suffix_splicing, ".RDS"))
      results_splicing <- tibble::as_tibble(results_splicing)
      results_expression <- tibble::as_tibble(results_expression, rownames = "gene_tx") %>%
        tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") 
    }
    
    if(!new & !limma){
      colnames(results_expression) <- colnames(results_expression) %>%
        stringr::str_replace("log2FoldChange", "lfc_deseq") %>%
        stringr::str_replace_all("deseq", "expression")
      colnames(results_splicing) <- colnames(results_splicing) %>%
        stringr::str_replace("log2fold_C_B", "lfc_dexseq") %>%
        stringr::str_replace("lfcSE", "lfc_dexseq") %>%
        stringr::str_replace_all("dexseq", "splicing") %>%
        stringr::str_replace_all("groupID", "gene") %>%
        stringr::str_replace_all("featureID", "transcript")
    }
    
    if("featureID" %in% colnames(results_splicing)) colnames(results_splicing) <- stringr::str_replace(colnames(results_splicing), "featureID", "transcript")
    
    
    genes_to_keep_splicing <- results_splicing %>%
      dplyr::group_by(gene) %>% 
      dplyr::summarise(n_significant = sum(padj < 0.05, na.rm = TRUE)) %>% 
      dplyr::filter(n_significant > 1) %>% 
      dplyr::pull(gene)
    genes_to_keep_expression <- results_expression %>%
      dplyr::group_by(gene) %>% 
      dplyr::summarise(n_significant = sum(padj < 0.05, na.rm = TRUE)) %>% 
      dplyr::filter(n_significant > 1) %>% 
      dplyr::pull(gene)
    
    genes_to_keep <- union(genes_to_keep_splicing, genes_to_keep_expression)
    
    results_joined <- results_expression %>% 
      dplyr::left_join(results_splicing, by = c("gene", "transcript"), suffix = c("_expression", "_splicing")) %>% 
      dplyr::filter(gene %in% genes_to_keep) %>% 
      dplyr::mutate(experiment = experiment_title) %>% 
      dplyr::select(-dplyr::starts_with("count"))
    
    
    concatenated_results <- concatenated_results %>% 
      dplyr::bind_rows(results_joined)
  }
  saveRDS(concatenated_results, paste0("results/concatenated", suffix,".RDS"))
  return(concatenated_results)
}
# concatenated_results <- concatenate_results(experiments)

# Concat gene level ----
concatenate_genes <- function(experiments, suffix = "_genes", new = FALSE, limma = FALSE){
  
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
    
    if(!new & !limma){
      colnames(aggregated_pvals) <- colnames(aggregated_pvals) %>%
        stringr::str_replace_all("deseq", "expression") %>%
        stringr::str_replace_all("dexseq", "splicing")
      aggregated_pvals <- aggregated_pvals %>% 
        dplyr::mutate(padj_expression = p.adjust(pvalue_expression, "fdr"),
                      padj_splicing = p.adjust(pvalue_splicing, "fdr"))
    }
    
    aggregated_pvals$experiment <- experiment_title
    
    
    concatenated_genes <- concatenated_genes %>% 
      dplyr::bind_rows(aggregated_pvals)
  }
  
  saveRDS(concatenated_genes, paste0("results/concatenated", suffix,".RDS"))
  return(concatenated_genes)
}
# concatenated_genes <- concatenate_genes(experiments)

concatenate_no_sva_genes <- function(experiments, suffix = "_no_sva_genes"){
  
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
  
  saveRDS(concatenated_no_sva_genes, paste0("results/concatenated", suffix,".RDS"))
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