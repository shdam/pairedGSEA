## Collection of plots for paper

# Initiatlize  ----
if(FALSE){
  pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
  library(tidyverse)
  theme_set(theme_classic(base_size = 20))

                    # Both  DGE          DGS        Neither
  color_scheme <- c("gray", "lightblue", "darkred", "black")
  fill_scale <- function(reorder = NULL) {
    if(!is.null(reorder)) color_scheme <- color_scheme[reorder]
    scale_fill_manual(values = color_scheme, 
                      na.value = "white")}
  color_scale <- function(reorder = NULL) {
    if(!is.null(reorder)) color_scheme <- color_scheme[reorder]
    scale_color_manual(values = color_scheme, 
                       na.value = "white")}
  
  
  utils::globalVariables(add = TRUE, c("color_scheme"))
}
# Figure 1 ----

## Figure 1B: Number of genes falsely significant when not including SVAs


plot_false_sva <- function(concatenated_no_sva_genes, concatenated_sva_genes){
  
  # Distinguish up and dowm-regulation
  gene_change_sva <- concatenated_sva_genes %>% 
    mutate(gene_change = paste0(sign(lfc_deseq))) %>% 
    dplyr::select(gene_change, experiment, gene)
  gene_change_no_sva <- concatenated_no_sva_genes %>% 
    mutate(gene_change = paste0(sign(lfc))) %>% 
    dplyr::select(gene_change, experiment, gene)
  
  # Number of genes falsely significant when not including SVAs
  false_sig <- gene_change_no_sva %>% 
    anti_join(gene_change_sva, by = c("experiment", "gene")) %>% 
    dplyr::count(experiment)
  
  
  plt_false_sig <- false_sig %>% 
    ggplot(aes(x = n)) +
    geom_density(alpha = 0.5, fill = "gray") +
    geom_vline(xintercept = median(false_sig$n, na.rm = TRUE), linetype = "dashed") +
    labs(
      x = paste0("Number of genes associated with surrogate variables\n", "(Median: ", round(median(false_sig$n, na.rm = TRUE), 3), ")"),
      y = "Density") +
    scale_x_log10() +
    theme(legend.position = "none")
  
  return(plt_false_sig)
  
}

if(FALSE){
  
  concatenated_no_sva_genes <- readRDS("results/concatenated_no_sva_genes.RDS")
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  concatenated_sva_genes <- concatenated_genes %>% filter(padj_deseq < 0.05)
  
  false_sva <- plot_false_sva(concatenated_no_sva_genes, concatenated_sva_genes)
  
  ggsave("figs/false_sva.png", false_sva)
  
}

# Figure 1C: Genes falsely significant as a fraction of # significant when not including SVAs
# Aka actual false discovery rates when not including SVAs

plot_fdr <- function(concatenated_no_sva_genes, concatenated_sva_genes){
  # Compute false discovery rate
  fdr <- gene_change_no_sva %>% 
    # anti_join(gene_change_sva, by = c("experiment", "gene")) %>% 
    dplyr::count(experiment) %>% 
    left_join(gene_change_no_sva %>% 
                inner_join(gene_change_sva, by = c("experiment", "gene")) %>% 
                dplyr::count(experiment), by = "experiment", suffix = c("_no_sva", "_overlap")) %>% 
    mutate(fdr = (n_no_sva - n_overlap) / n_no_sva)
  
  plt_fdr <- fdr %>% 
    ggplot(aes(x = fdr)) +
    geom_histogram(fill = "darkgray", color = "white", alpha = 0.8) +
    geom_vline(xintercept = median(fdr$fdr, na.rm = TRUE), linetype = "dashed") +
    labs(x = paste0("Expected False Discovery Rate\nWithout Proper Correction for Confounders", "\n(Median: ", round(median(fdr$fdr, na.rm = TRUE), 3), ")"),
         y = "Count") +
    theme(legend.position = "none")
  
  return(plt_fdr)
}

if(FALSE){
  # Load data
  concatenated_no_sva_genes <- readRDS("results/concatenated_no_sva_genes.RDS")
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  concatenated_sva_genes <- concatenated_genes %>% filter(padj_deseq < 0.05)
  
  
  fig1c <- plot_fdr(concatenated_no_sva_genes, concatenated_sva_genes)
  ggsave("figs/fig1c.png", fig1c)
  
}

# Figure 2 ----
# Figure 2A: Genes per experiment

plot_gene_counts <- function(concatenated_genes){
  # Count number of significant genes
  found_genes <- concatenated_genes %>% 
    group_by(experiment) %>%
    summarise(`Differential\nExpression` = sum(padj_deseq < 0.05, na.rm = TRUE),
              `Differential\nSplicing` = sum(padj_dexseq < 0.05, na.rm = TRUE),
              Overlap = sum(padj_dexseq < 0.05 & padj_deseq < 0.05, na.rm = TRUE),
              `Only\nDifferential\nSplicing` = sum((padj_dexseq < 0.05) & !(padj_deseq < 0.05), na.rm = TRUE))
  ## Plot as sina
  plt_gene_counts <- found_genes %>% 
    pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Genes") %>% 
    mutate(Analysis = factor(Analysis, levels = c("Differential\nExpression",
                                                  "Differential\nSplicing",
                                                  "Overlap",
                                                  "Only\nDifferential\nSplicing"))) %>% 
    ggplot(aes(y = Genes, x = Analysis, color = Analysis)) +
    scale_y_log10() +
    ggforce::geom_sina() +
    geom_boxplot(width = 0.05) +
    color_scale(reorder = c(2, 3, 1, 4)) +
    labs(x = "",
         y = "Significant genes\nper experiment") +
    theme(legend.position = "none")
  return(plt_gene_counts)
  
}


if(FALSE){
  # Load data
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  
  plt_gene_counts <- plot_gene_counts(concatenated_genes)
  ggsave(plot = plt_gene_counts, filename = "figs/gene_counts.png")
  
}


# Figure 2B: Fraction of significant genes

plot_gene_fraction <- function(concatenated_genes){
  
  # Compute gene fractions
  gene_fractions <- concatenated_genes %>% 
    group_by(experiment) %>% 
    summarise(deseq = sum(padj_deseq < 0.05, na.rm = TRUE),
              num_genes = n(),
              dexseq = sum(padj_dexseq < 0.05, na.rm = TRUE)) %>% 
    mutate(`Differential\nExpression` = deseq/num_genes,
            `Differential\nSplicing` = dexseq/num_genes) %>% 
    pivot_longer(cols = starts_with("Differential"), names_to = "Method", values_to = "Fraction")
  
  # Sina plot of fraction of found genes
  plt_gene_fractions <-  gene_fractions %>% 
    ggplot() +
    aes(x = Method, y = Fraction, color = Method) +
    ggforce::geom_sina() +
    geom_boxplot(width = 0.05) +
    color_scale(reorder = c(2, 3)) +
    labs(x = "",
         y = "Fraction of significant genes") +
      theme(legend.position = "none")
  
  return(plt_gene_fractions)
  
}

if(FALSE){
  # Load data
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  
  plt_gene_fractions <- plot_gene_fraction(concatenated_genes)
  ggsave(plot = plt_gene_fractions, filename = "figs/gene_fractions.png")
  
  
  # Number and fraction of significant genes
  concatenated_genes %>% 
    group_by(experiment) %>% 
    summarise(deseq = sum(padj_deseq < 0.05, na.rm = TRUE),
              num_genes = n(),
              dexseq = sum(padj_dexseq < 0.05, na.rm = TRUE)) %>% 
    summarise(med_deseq = median(deseq), med_dexseq = median(dexseq),
              frac_deseq = median(deseq/num_genes),
              frac_dexseq = median(dexseq/num_genes))
  
}


# Figure 2C + 2D: DGS and DGE Similarities

# Compute gene similarity
compute_gene_similarity <- function(concatenated_genes){
  gene_similarity <- concatenated_genes %>%
    mutate(DGE = padj_deseq < 0.05,
           DGS = padj_dexseq < 0.05) %>% 
    group_by(experiment) %>% 
    summarise(Similarity = proxy::simil(x = DGE, y = DGS, 
                                        by_rows = FALSE, method = "Simpson")[1],
              n_DGE = sum(DGE, na.rm = T),
              n_DGS = sum(DGS, na.rm = T),
              p_overlap = (min(n_DGS, n_DGE)*Similarity) / n_DGE,
              p_signal = n_DGS / (n_DGE + n_DGS - min(n_DGS, n_DGE)*Similarity)
    ) %>% 
    filter(!is.na(Similarity))
  
  return(gene_similarity)
}

plot_dgs_affected_dge <- function(gene_similarity){
  
  median_overlap <- median(gene_similarity$p_overlap, na.rm = TRUE)
  
  plt_affected <- gene_similarity %>% 
    ggplot(aes(x = p_overlap)) +
    geom_density(alpha = 0.5, fill = "darkgray") +
    geom_vline(xintercept = median_overlap, linetype = "dashed") +
    labs(
      x = paste0("Differentially expressed genes\naffected by differential splicing",
                 "\n(Median: ", round(median_overlap, 3), ")"),
      y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    theme(legend.position = "none")
  
  return(plt_affected)
}

plot_dgs_signal <- function(gene_similarity){
  
  median_signal <- median(gene_similarity$p_signal, na.rm = TRUE)
  
  plt_dgs_signal <- gene_similarity %>% 
    ggplot(aes(x = p_signal)) +
    geom_density(alpha = 0.5, fill = "darkgray") +
    geom_vline(xintercept = median_signal, linetype = "dashed") +
    labs( x = paste0("Fraction of total signal mediated by differential splicing",
                     "\n(Median: ", round(median_signal, 3), ")"),
          y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    theme(legend.position = "none")
  
  return(plt_dgs_signal)
}


  
if(FALSE){
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  
  gene_similarity <- compute_gene_similarity(concatenated_genes)
  
  # Medians
  median(gene_similarity$p_overlap)
  median(gene_similarity$p_signal)
  
  dgs_affected <- plot_dgs_affected_dge(gene_similarity)
  dgs_signal <- plot_dgs_signal(gene_similarity)
  
  ggsave("figs/dgs_affected.png", dgs_affected)
  ggsave("figs/dgs_signal.png", dgs_signal)
  
}


# Add TPM information ----
if(FALSE){
  
  source("tmp/run_experiment.R")
  tpm_file = "/home/databases/archs4/v11/human_transcript_v11_tpm.h5"
  gtf <- readRDS("gtfextract.rds")
  gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)
  
  extract_tpm <- function(exp, both_genes){
    print(exp)
    experiment_title <- str_remove(exp, "\\d+_GSE\\d+_") %>% 
      str_remove("GSE\\d+_")
    dataset <- str_remove(exp, fixed(paste0("_", experiment_title)))
    
    md <- readxl::read_xlsx(paste0("metadata/", dataset, ".xlsx"))
    comparison <- md$`comparison (baseline_v_condition)`[which(md$`comparison_title (empty_if_not_okay)` == experiment_title)] %>% str_split("v", simplify = TRUE)
    
    md <- md %>% 
      filter(group_nr %in% comparison)
    
    keep_genes <- filter(both_genes, experiment == exp) %>% 
      pull(gene)
    
    
    tx_tpm <- load_archs4(md$id, tpm_file, gtf = gtf)
    
    colnames(tx_tpm) <- paste0(md$id, "_", md$group_nr)
    
    tx_tpm <- tx_tpm %>% 
      as_tibble(rownames = "gene_tx") %>% 
      separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
      filter(gene %in% keep_genes) %>%
      mutate(experiment = exp,
             tpm_baseline = rowMeans(select(., ends_with(comparison[1]))),
             tpm_condition = rowMeans(select(., ends_with(comparison[2])))) %>% 
      select(-starts_with("GSM"))
    
    return(tx_tpm)
  }
  
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  both_genes <- concatenated_genes %>% 
    filter(padj_dexseq < 0.05 & padj_deseq < 0.05,
           gene != "NA")
  
  tpms <- lapply(unique(both_genes$experiment), extract_tpm, both_genes) %>% 
    bind_rows() %>% 
    as_tibble()
  
  saveRDS(tpms, file = "results/tpms.RDS")
  
  
}

# Prepare TPM plots

if(FALSE){
  
  tpms <- readRDS("results/tpms.RDS")
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  both_genes <- concatenated_genes %>% 
    filter(padj_dexseq < 0.05 & padj_deseq < 0.05,
           gene != "NA")
  
  dgs_isoforms <- readRDS("results/concatenated_results.RDS") %>% 
    semi_join(both_genes, by = c("gene", "experiment"))
  
}

plot_switch_fraction <- function(dgs_isoforms) {
  
  count_switches <- dgs_isoforms %>% 
    filter(padj_dexseq < 0.05) %>% 
    mutate(sign = sign(log2FC_dexseq)) %>% 
    distinct(sign, gene, experiment) %>% 
    count(gene, experiment) %>% 
    filter(n > 1) %>% 
    count(experiment, name = "events") %>% 
    left_join(both_genes %>% count(experiment, name = "num_genes"), by = c("experiment")) %>% 
    mutate(percent = events/num_genes)
  
  
  
  median_switches <- median(count_switches$percent, na.rm = TRUE)
  
  plt_switches <- count_switches %>% 
    ggplot() +
    aes(x = percent) +
    geom_density(alpha = 0.5, fill = "darkgray") +
    geom_vline(xintercept = median_switches, linetype = "dashed") +
    labs(x = paste0("Fraction of genes with an isoform switch",
                    "\n(Median: ", round(median_switches, 3), ")"),
         y = "Density") +
    coord_cartesian(xlim = c(0, 1))
  
  return(plt_switches)
}


if(FALSE){
  
  
  plt_switches <- plot_switch_fraction(dgs_isoforms)
  plt_switches
  ggsave("figs/isoform_switches.png", plt_switches)
}

count_switched_majority <- function(dgs_isoforms, tpms){
  # Define LFC direction and add tpm information
  switched_majorities <- dgs_isoforms %>% 
    filter(padj_dexseq < 0.05) %>% 
    mutate(sign = as.character(sign(log2FC_dexseq))) %>% 
    left_join(tpms, by = c("gene", "transcript", "experiment")) %>% 
    group_by(gene, experiment) %>% 
    # Identify switches and transcript majorities
    summarise(major_baseline = transcript[which(tpm_baseline == max(tpm_baseline))], 
              major_condition = transcript[which(tpm_condition == max(tpm_condition))],
              switch = "-1" %in% sign & "1" %in% sign) %>% 
    filter(switch) %>% 
    group_by(experiment) %>% 
    # Count changes in majority
    summarise(num_genes = n(), major_change = sum((major_baseline != major_condition))) %>% 
    mutate(change_percent = major_change/num_genes)
  
  return(switched_majorities)
}

plot_switch_majority <- function(switched_majorities){
  
  median_switched_majority <- median(switched_majorities$change_percent, na.rm = TRUE)
  
  plt_switched_majorities <- switched_majorities %>%
    ggplot() +
    aes(x = change_percent) +
    geom_density(alpha = 0.5, fill = "darkgray") +
    geom_vline(xintercept = median_switched_majority, linetype = "dashed") +
    labs(x = paste0("Fraction of isoform switches leading to\nchange in most expressed isoform",
                    "\n(Median: ", round(median_switched_majority, 3), ")"),
         y = "Density") +
    coord_cartesian(xlim = c(0, 1))
  
  
  return(plt_switched_majorities)
  
}


if(FALSE){
  
  switched_majorities <- count_switched_majority(dgs_isoforms, tpms)
  
  plt_switch_majority <- plot_switch_majority(switched_majorities)
  plt_switch_majority
  ggsave("figs/switch_majority.png", plt_switch_majority)
}


# Figure 3 ----

# Figure 3A: Significant pathways per experiment

plot_pathway_count <- function(ora_all){
  
  found_pathways <- ora_all %>% 
    group_by(experiment) %>% 
    summarise(`Differential\nExpression` = sum(padj_deseq < 0.05, na.rm = TRUE),
              `Differential\nSplicing` = sum(padj_dexseq < 0.05, na.rm = TRUE),
              Overlap = sum((padj_dexseq < 0.05) & (padj_deseq < 0.05), na.rm = TRUE),
              `Only\nDifferential\nSplicing` = sum(!(padj_deseq < 0.05) & (padj_dexseq < 0.05), na.rm = TRUE)
    )
  
 # Sina plot of found pathways
  plt_pathway_counts <- found_pathways %>% 
    pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Pathways") %>% 
    mutate(Analysis = factor(Analysis, levels = c("Differential\nExpression",
                                                  "Differential\nSplicing",
                                                  "Overlap",
                                                  "Only\nDifferential\nSplicing")
    )) %>% 
    ggplot() +
    aes(y = Pathways, x = Analysis, color = Analysis) +
    ggforce::geom_sina() +
    scale_y_log10() +
    geom_boxplot(width = 0.05) +
    color_scale(reorder = c(2, 3, 1, 4)) + 
    labs(x = "",
         y = "Significant pathways\nper experiment") +
    theme(legend.position = "none")
  
  return(plt_pathway_counts)
  
}

if(FALSE){
  
  ora_all <- readRDS("results/ora_all.RDS")
  
  plt_pathway_count <- plot_pathway_count(ora_all)
  plt_pathway_count
  ggsave(plot = plt_pathway_count, filename = "figs/pathway_count.png")
}

# Figure 3B: Telomere pathways ----

if(FALSE){
  ora_all <- readRDS("results/ora_all.RDS")
  exp <- "77_GSE61220_TNF Treatment 12hrs"
  ora <- ora_all %>% 
    filter(experiment == exp)
  
  
  telomere_pathways <- plot_ora(ora, pattern = "Telomer")
  ggsave(plot = telomere_pathways, filename = "figs/telomere_pathways.png")
}

# Figure 3C: Repair pathways----

if(FALSE){
  ora_all <- readRDS("results/ora_all.RDS")
  exp <- "79_GSE139262_SMARCB1 overexpression"
  ora <- ora_all %>% 
    filter(experiment == exp)
  
  
  repair_pathways <- plot_ora(ora, pattern = "Repair")
  ggsave(plot = repair_pathways, filename = "figs/repair_pathways.png")
}

# Figure 3D: A figure summarizing correlations off gene-set enrichments across all datasets

plot_ora_correlation <- function(ora_all, threshold = 0){
  
  correlation <- ora_all %>% 
    filter(padj_dexseq < 0.05 | padj_deseq < 0.05,
           !is.na(padj_dexseq),
           !is.na(padj_deseq)) %>% 
    group_by(experiment) %>% 
    anti_join(count(.) %>% filter(n < threshold), by = "experiment") %>%
    summarise(correlation = cor(enrichment_score_dexseq,
                                enrichment_score_deseq, method = "spearman")) %>% 
    filter(!is.na(correlation))
  
  median_correlation <- median(correlation$correlation)
  
  plt_correlation <- correlation %>% 
    ggplot(aes(x = correlation)) +
    geom_histogram(fill = "darkgray", color = "white", alpha = 0.8) +
    geom_vline(xintercept = median_correlation, linetype = "dashed") +
    labs(x = paste0("Spearman's ρ between gene-set enrichment scores", "\n(Median: ", round(median_correlation, 3), ")"),
         y = "Count")
   return(plt_correlation)
}

plot_ora_correlation_facet <- function(ora_all, threshold = 50){
  correlation <- ora_all %>% 
    filter(!is.na(padj_dexseq),
           !is.na(padj_deseq)) %>% 
    mutate(association = case_when(padj_dexseq < 0.05 & padj_deseq < 0.05 ~ "Both",
                                   padj_dexseq < 0.05 ~ "DGS",
                                   padj_deseq < 0.05 ~ "DGE",
                                   TRUE ~ "Neither"),
           association = factor(association, levels = c("Both", "DGE", "DGS", "Neither"))) %>% 
    filter(association != "Neither") %>% 
    group_by(experiment, association) %>% 
    anti_join(count(.) %>% filter(n < threshold), by = c("experiment", "association")) %>%
    summarise(correlation = cor(enrichment_score_dexseq,
                                enrichment_score_deseq, method = "spearman")) %>% 
    filter(!is.na(correlation), !is.na(association))
  
  # Medians
  medians <- correlation %>% 
    group_by(association) %>% 
    summarise(median = median(correlation, na.rm = TRUE))
  print(medians)
  
  plt_correlation <- correlation %>% 
    # filter(association != "Neither") %>% 
    ggplot(aes(x = correlation, fill = association)) +
    geom_histogram(color = "white") +
    geom_vline(
      data = . %>%
        group_by(association) %>%
        summarise(line = median(correlation)),
      mapping = aes(xintercept = line, color = association),
      linetype = "dashed"
    ) +
    labs(x = "Spearman's ρ between gene-set enrichment scores",
         y = "Count") +
    facet_grid(~association, space="free_x")  +
    fill_scale() +
    color_scale() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle=-90))

  return(plt_correlation)
}

if(FALSE){
  
  ora_all <- readRDS("results/ora_all.RDS")
  
  ora_correlation <- plot_ora_correlation(ora_all)
  ora_correlation
  ggsave(plot = ora_correlation, filename = "figs/ora_correlation.png")
  ora_correlation_facet <- plot_ora_correlation_facet(ora_all)
  ora_correlation_facet
  ggsave(plot = ora_correlation_facet, filename = "figs/ora_correlation_facet.png")
}


# Figure 3E: A figure summarizing relative risk changes across all datasets ----


compute_rr_shifts <- function(ora_all){
  
  rr_shifts <- ora_all %>% 
    dplyr::mutate(Shift = relative_risk_dexseq - relative_risk_deseq,
                  association = dplyr::case_when(
                    padj_dexseq < 0.05 & padj_deseq < 0.05 ~ "Both",
                    padj_dexseq < 0.05 ~ "DGS",
                    padj_deseq < 0.05 ~ "DGE",
                    TRUE ~ "Neither"),
                  association = factor(association, levels = c("Both", "DGE", "DGS", "Neither"))) %>% 
    dplyr::filter(association != "Neither") %>%
    dplyr::arrange(desc(Shift))
  
  return(rr_shifts)
}



plot_rr_ridges <- function(rr_shifts){
  
  rr_ridges <- rr_shifts %>%
    group_by(experiment, association) %>% 
    summarise(Shift = median( abs(Shift) / (min(relative_risk_dexseq, relative_risk_deseq))) * 100) %>%
    filter(!is.infinite(Shift),
           !is.na(Shift)) 
  
  rr_ridges_medians <- rr_ridges %>% 
    group_by(association) %>% 
    summarise(median = median(Shift))
  print(rr_ridges_medians)
  
  plt_rr_ridges <- rr_ridges %>% 
    ggplot(aes(x = Shift, y = association, color = association, fill = association)) +
    ggridges::geom_density_ridges(scale = 4, alpha = 0.4, quantile_lines = TRUE, quantiles = 2) +
    color_scale() +
    fill_scale() +
    labs(x = "Shift in gene-set relative risk as a percentage\nof the smaller risk",
         y = ""
    ) +
    scale_x_log10() +
    theme(legend.position = "none")
  
  return(plt_rr_ridges)
}

plot_rr_median <- function(rr_shifts){
  
  rr_shift_median <- rr_shifts %>%
    group_by(experiment) %>% 
    summarise(Shift = median( abs(Shift) / min( relative_risk_dexseq, relative_risk_deseq) ) * 100) %>% 
    filter(!is.infinite(Shift), !is.na(Shift))
  
  median_of_median <- median(rr_shift_median$Shift)
  plt_rr_median <- rr_shift_median %>% 
    ggplot(aes(x = Shift)) +
    geom_histogram(fill = "darkgray", color = "white", alpha = 0.8) +
    geom_vline(xintercept = median_of_median, linetype = "dashed") +
    labs(x = paste0("Shift in gene-set relative risk"," as a percentage\nof the smaller risk", " (Median: ", round(median_of_median, 2), ")"),
         y = "Count"
    ) +
    scale_x_log10() +
    theme(legend.position = "none")
  
  return(plt_rr_median)
}


if(FALSE){
  
  ora_all <- readRDS("results/ora_all.RDS")
  rr_shifts <- compute_rr_shifts(ora_all)
  
  
  rr_ridges <- plot_rr_ridges(rr_shifts)
  rr_ridges
  ggsave(plot = rr_ridges, filename = "figs/rr_ridges.png")
  rr_median <- plot_rr_median(rr_shifts)
  rr_median
  ggsave(plot = rr_median, filename = "figs/rr_median.png")
}

