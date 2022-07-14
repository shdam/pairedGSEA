
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
    
    genes_to_keep_dexseq <- results_dexseq %>% 
      dplyr::group_by(groupID) %>% 
      dplyr::summarise(n_significant = sum(padj < 0.05, na.rm = TRUE)) %>% 
      dplyr::filter(n_significant > 1) %>% 
      dplyr::pull(groupID)
    genes_to_keep_deseq <- results_deseq %>% 
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


# Concat SVA level ----
concatenate_design <- function(experiments){
  
  concatenaded_design <- tibble::tibble()
  
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
    
    
    concatenaded_design <- concatenaded_design %>% 
      dplyr::bind_rows(des)
  }
  saveRDS(concatenaded_design, "results/concatenaded_design.RDS")
  return(concatenaded_design)
}

concatenate_sva_genes <- function(experiments){
  
  concatenated_sva_genes <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    row <- experiments[num, ]
    ### Load metadata
    md_file <- row$filename
    data_name <- basename(md_file) %>% 
      stringr::str_remove(".xlsx") %>% 
      stringr::str_remove("csv") 
    ### Define experiment details
    baseline_case <- row$`comparison (baseline_v_condition)` %>% stringr::str_split(pattern = "v", simplify = TRUE) %>% as.character()
    experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
    ### Check that results exists
    message("Adding ", experiment_title)
    dds <- readRDS(paste0("results/", experiment_title, "_dds.RDS"))
    
    svs <- as.character(DESeq2::design(dds))[2] %>% 
      stringr::str_count("sv")
    
    found_genes <- c()
    
    for(sv in 1:length(svs)){
      
      deseq_results <- DESeq2::results(dds,
                                       name = paste0("sv", sv),
                                       parallel = FALSE)
      
      # Convert result to tibble
      deseq_results <- deseq_results %>% 
        tibble::as_tibble(rownames = "gene_tx") %>% 
        tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
        dplyr::rename(log2FC_deseq = log2FoldChange)
      
      # Aggregate to gene level
      deseq_aggregated <- aggregate_pvalue(deseq_results, gene = "gene", weights = "baseMean", lfc = "log2FC_deseq", type = "deseq")
      
      # Extract found genes
      found_genes <- found_genes %>% 
        union(
          deseq_aggregated %>% 
          dplyr::mutate(padj = stats::p.adjust(pvalue, "fdr")) %>% 
          dplyr::filter(padj < 0.05) %>% 
          dplyr::pull(gene)
        )
      
    }
    
    
    concatenated_sva_genes <- concatenated_sva_genes %>% 
      dplyr::bind_rows(tibble("gene" = found_genes, "experiment" = experiment_title))
  }
  saveRDS(concatenated_sva_genes, "results/concatenated_sva_genes.RDS")
  return(concatenated_sva_genes)
}


run_experiment_no_sva <- function(row, archs4db = NULL, tx_count = NULL, group_col = "group_nr",
                                  tpm = TRUE, prefilter = 10, parallel = TRUE, gtf = NULL, quiet = FALSE){
  
  if(typeof(row) == "character"){ # Convert apply-made row to tibble
    row <- tibble::as_tibble(row, rownames = "names") %>% 
      tidyr::pivot_wider(values_from = value, names_from = names)
  }
  
  
  message("Running on ", row$study)
  
  ### Load metadata
  md_file <- row$filename
  data_name <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove(".csv") 
  
  ### Define tpm file
  if(tpm) tpm <- stringr::str_replace(archs4db, "counts", "tpm")
  ### Define experiment details
  baseline_case <- row$`comparison (baseline_v_condition)` %>% stringr::str_split(pattern = "v", simplify = TRUE) %>% as.character()
  experiment_title <- paste0(data_name, "_", row$`comparison_title (empty_if_not_okay)`)
  
  
  ### Prepare for DE
  if(is.null(tx_count)){
    tx_count <- prepare_tx_count(
      metadata = md_file,
      gtf = gtf,
      archs4db = archs4db,
      group_col = group_col,
      baseline_case = baseline_case
    )
  }
  
  baseline <- baseline_case[1]
  case <- baseline_case[2]
  sample_col <- "id"
  
  metadata <- prepare_metadata(md_file, group_col, paste(c(baseline, case)))
  metadata <- metadata[metadata[[sample_col]] %in% colnames(tx_count), ]
  
  tx_count <- tx_count[, metadata[[sample_col]]]
  
  ### Prefiltering
  if(prefilter) tx_count <- pre_filter(tx_count, prefilter)
  
  if(!quiet) message("Converting count matrix to DESeqDataSet")
  # Create DDS from count matrix
  design <- formularise_vector(c(group_col))
  dds <- convert_matrix_to_dds(tx_count, metadata, design)
  
  
  deseq_results <- run_deseq(
    dds,
    group_col = group_col,
    baseline = baseline,
    case = case,
    experiment_title = experiment_title
  )
  
  deseq_aggregated <- aggregate_pvalue(deseq_results, gene = "gene", weights = "baseMean", lfc = "log2FC_deseq", type = "deseq") %>% 
    dplyr::mutate(padj = stats::p.adjust(pvalue, "fdr"),
                  experiment = experiment_title) %>% 
    dplyr::filter(padj < 0.05)
  
  
  
  return(deseq_aggregated)
  
}

deseq_no_sva <- function(experiments, archs4db, gtf){
  
  concatenated_no_sva_genes <- tibble::tibble()
  
  for(num in 1:nrow(experiments)){
    deseq_genes <- run_experiment_no_sva(experiments[num, ], archs4db, gtf = gtf)
    
    concatenated_no_sva_genes <- concatenated_no_sva_genes %>% 
      bind_rows(deseq_genes)
  }
  saveRDS(concatenated_no_sva_genes, "results/concatenated_no_sva_genes.RDS")
  return(concatenated_no_sva_genes)
}

# deseq_genes_no_sva <- deseq_no_sva(experiments, archs4db, gtf)


concatenate_sva_correlation <- function(experiments){
  
  concatenaded_correlation <- tibble::tibble()
  
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
    
    concatenaded_correlation <- concatenaded_correlation %>% 
      dplyr::bind_rows(correlation)
  }
  saveRDS(concatenaded_correlation, "results/concatenate_sva_correlation.RDS")
  return(concatenaded_correlation)
}
# cors <- concatenaded_correlation %>% pivot_longer(cols = starts_with("sv"), names_to = "sv", values_to = "correlation")