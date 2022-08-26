## Collection of plots for paper

initialize <- TRUE
fig1 <- TRUE
fig1_init <- TRUE
fig2 <- TRUE
fig2_init <- TRUE
fig3 <- TRUE
fig_format <- "png"

# Initialize  ----
if(initialize){
  pkgload::load_all(path = "/home/projects/shd/pairedGSEA")
  library(tidyverse)
  library(patchwork)
  theme_set(theme_classic(base_size = 19))

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
  
  
  full_figure_theme <- function(){
      theme(plot.tag = element_text(face = "plain", size = 10))
  }
  
  utils::globalVariables(add = TRUE, c("color_scheme"))
  
  
  # Save metadata overview
  
  source("tmp/run_experiment.R")
  experiments <- combine_experiments()
  experiments %>% select(-c("final_description", "group_nr", "id", "filename")) %>% 
    writexl::write_xlsx("figs/Experiments.xlsx")

  
}
# Figure 1 ----

## Figure 1B: Number of genes falsely significant when not including SVAs


plot_false_sva <- function(gene_change_no_sva, gene_change_sva){
  
  
  
  # Number of genes falsely significant when not including SVAs
  false_sig <- gene_change_no_sva %>% 
    anti_join(gene_change_sva, by = c("experiment", "gene")) %>%
    dplyr::count(experiment)
  
  
  plt_false_sig <- false_sig %>% 
    ggplot(aes(x = n)) +
    geom_density(alpha = 0.5, fill = "gray") +
    geom_vline(xintercept = median(false_sig$n, na.rm = TRUE), linetype = "dashed") +
    labs(
      x = paste0("Number of genes associated with confounders\n", "(Median: ", round(median(false_sig$n, na.rm = TRUE), 3), ")"),
      y = "Density") +
    scale_x_log10() +
    theme(legend.position = "none")
  
  return(plt_false_sig)
  
}

# Figure 1C: Genes falsely significant as a fraction of # significant when not including SVAs
# Aka actual false discovery rates when not including SVAs

plot_fdr <- function(gene_change_no_sva, gene_change_sva){
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
    labs(x = paste0("Expected false discovery rate\nwithout proper correction for confounders", "\n(Median: ", round(median(fdr$fdr, na.rm = TRUE), 3), ")"),
         y = "Count") +
    theme(legend.position = "none")
  
  return(plt_fdr)
}


if(fig1_init){
  
  concatenated_no_sva_genes <- readRDS("results/concatenated_no_sva_genes.RDS")
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  concatenated_sva_genes <- concatenated_genes %>% filter(padj_deseq < 0.05)
  
  concatenated_no_sva_genes <- concatenated_no_sva_genes %>% 
    dplyr::mutate(experiment = stringr::str_remove(experiment, "_no_sva"))
  
  # Distinguish up and dowm-regulation
  gene_change_sva <- concatenated_sva_genes %>% 
    dplyr::filter(padj_deseq < 0.05) %>% 
    mutate(gene_change = paste0(sign(lfc_deseq))) %>% 
    dplyr::select(gene_change, experiment, gene)
  gene_change_no_sva <- concatenated_no_sva_genes %>% 
    dplyr::filter(padj < 0.05) %>% 
    mutate(gene_change = paste0(sign(lfc))) %>% 
    dplyr::select(gene_change, experiment, gene)
  
}

if(fig1){
  
  false_sva <- plot_false_sva(gene_change_no_sva, gene_change_sva)
  
  ggsave(paste0("figs/false_sva.", fig_format), false_sva)
  
  fig1c <- plot_fdr(gene_change_no_sva, gene_change_sva)
  ggsave(paste0("figs/fig1c.", fig_format), fig1c)
  
  ## Full figure ----
  
  
  theme_set(theme_classic(base_size = 8))
  # magick
  library(cowplot)
  library(magick)
  
  img_g <- ggdraw() +
    draw_image("figs/overview.png") +
    draw_plot(ggplot() + theme_void())
  
  
  layout <- "
  AAAA
  AAAA
  BBCC
  "
  Figure1 <- img_g + false_sva + fig1c +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') & 
    full_figure_theme()
  
  # Figure1
  
  ggsave(paste0("figs/Figure1.", fig_format), Figure1, width = 6, height = 2)
  theme_set(theme_classic(base_size = 19))
}







# Figure 2 ----
# Figure 2A: Genes per experiment

plot_gene_counts <- function(concatenated_genes){
  # Count number of significant genes
  found_genes <- concatenated_genes %>% 
    group_by(experiment) %>%
    summarise(DGE = sum(padj_deseq < 0.05, na.rm = TRUE),
              DGS = sum(padj_dexseq < 0.05, na.rm = TRUE),
              Overlap = sum(padj_dexseq < 0.05 & padj_deseq < 0.05, na.rm = TRUE),
              Only_DGS = sum((padj_dexseq < 0.05) & !(padj_deseq < 0.05), na.rm = TRUE))
  ## Plot as sina
  plt_gene_counts <- found_genes %>% 
    pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Genes") %>% 
    mutate(Analysis = factor(Analysis, levels = c("DGE",
                                                  "DGS",
                                                  "Overlap",
                                                  "Only_DGS")),
           Analysis =  recode(Analysis,
                              DGE = paste0("Differential\nExpression\n", "(Median: ", median(found_genes$DGE), ")"),
                              DGS = paste0("Differential\nSplicing\n", "(Median: ", median(found_genes$DGS), ")"),
                              Overlap = paste0("Overlap\n", "(Median: ", median(found_genes$Overlap), ")"),
                              Only_DGS = paste0("Only\nDifferential\nSplicing\n", "(Median: ", median(found_genes$Only_DGS), ")"))
           ) %>% 
    ggplot(aes(y = Genes, x = Analysis, color = Analysis)) +
    scale_y_log10() +
    ggforce::geom_sina() +
    geom_boxplot(width = 0.05) +
    color_scale(reorder = c(2, 3, 1, 4)) +
    labs(x = "",
         y = "Significant genes\nper comparison") +
    theme(legend.position = "none")
  return(plt_gene_counts)
  
}


if(fig2){
  # Load data
  concatenated_genes <- readRDS(paste0("results/concatenated_genes.", fig_format))
  
  plt_gene_counts <- plot_gene_counts(concatenated_genes)
  ggsave(plot = plt_gene_counts, filename = paste0("figs/gene_counts.", fig_format))
  
}


# Figure 2B: Fraction of significant genes

plot_gene_fraction <- function(concatenated_genes){
  
  # Compute gene fractions
  gene_fractions <- concatenated_genes %>% 
    group_by(experiment) %>% 
    summarise(deseq = sum(padj_deseq < 0.05, na.rm = TRUE),
              num_genes = n(),
              dexseq = sum(padj_dexseq < 0.05, na.rm = TRUE)) %>% 
    mutate(`DGE` = deseq/num_genes,
            `DGS` = dexseq/num_genes) 
    
  
  # Sina plot of fraction of found genes
  plt_gene_fractions <-  gene_fractions %>% 
    pivot_longer(cols = starts_with("DG"), names_to = "Method", values_to = "Fraction") %>% 
    mutate(Method = recode(factor(Method), 
                           DGE = paste0("Differential\nExpression\n", "(Median: ", round(median(gene_fractions$DGE), 3), ")"),
                           DGS = paste0("Differential\nSplicing\n", "(Median: ", round(median(gene_fractions$DGS), 3), ")"))
           ) %>% 
    ggplot() +
    aes(x = Method, y = Fraction, color = Method) +
    ggforce::geom_sina() +
    geom_boxplot(width = 0.05) +
    color_scale(reorder = c(2, 3)) +
    labs(x = "",
         y = "Fraction of significant genes\nof genes tested") +
      theme(legend.position = "none")
  
  return(plt_gene_fractions)
  
}

if(fig2){
  # Load data
  # concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  
  plt_gene_fractions <- plot_gene_fraction(concatenated_genes)
  ggsave(plot = plt_gene_fractions, filename = paste0("figs/gene_fractions.", fig_format))
  
  
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
    labs( x = paste0("Fraction of total signal\nmediated by differential splicing",
                     "\n(Median: ", round(median_signal, 3), ")"),
          y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    theme(legend.position = "none")
  
  return(plt_dgs_signal)
}


  
if(fig2){
  # concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  
  gene_similarity <- compute_gene_similarity(concatenated_genes)
  
  # Medians
  median(gene_similarity$p_overlap)
  median(gene_similarity$p_signal)
  
  dgs_affected <- plot_dgs_affected_dge(gene_similarity)
  dgs_signal <- plot_dgs_signal(gene_similarity)
  
  # ggsave("figs/dgs_affected_signal.png", dgs_affected/dgs_signal)
  
  ggsave(paste0("figs/dgs_affected.", fig_format), dgs_affected)
  ggsave(paste0("figs/dgs_signal.", fig_format), dgs_signal)
  
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

if(fig2_init){
  
  tpms <- readRDS("results/tpms.RDS")
  concatenated_genes <- readRDS("results/concatenated_genes.RDS")
  both_genes <- concatenated_genes %>% 
    filter(padj_dexseq < 0.05 & padj_deseq < 0.05,
           gene != "NA")
  
  dgs_isoforms <- readRDS("results/concatenated_results.RDS") %>% 
    semi_join(both_genes, by = c("gene", "experiment")) %>% 
    dplyr::rename(log2FC_dexseq = log2fold_C_B,
                  log2FC_deseq = log2FoldChange)
  
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
    labs(x = paste0("Fraction of \"Both\" genes with\nan isoform switch",
                    "\n(Median: ", round(median_switches, 3), ")"),
         y = "Density") +
    coord_cartesian(xlim = c(0, 1)) +
    theme(legend.position = "none")
  
  return(plt_switches)
}

compute_IF <- function(tpms, dgs_isoforms, both_genes){
  tpm_IF <- tpms %>% 
    right_join(dgs_isoforms, by = c("transcript", "gene", "experiment")) %>% 
    semi_join(both_genes, by = "gene") %>% 
    group_by(experiment, gene) %>% 
    mutate(IF_baseline = tpm_baseline / sum(tpm_baseline),
           IF_condition = tpm_condition / sum(tpm_condition),
           dIF = IF_condition - IF_baseline)
  
  tpm_changing_IF <- tpm_IF %>% 
    filter(padj_dexseq < 0.05) %>% 
    summarise(changing_IF_baseline = sum(IF_baseline, na.rm = TRUE),
              changing_IF_condition = sum(IF_condition, na.rm = TRUE),
              changing_IF_mean = mean(c(changing_IF_baseline, changing_IF_condition))) %>% 
    group_by(experiment) %>% 
    summarise(median_IF_baseline = median(changing_IF_baseline),
              mean_baseline = mean(changing_IF_baseline),
              median_condition = median(changing_IF_condition),
              mean_condition = mean(changing_IF_condition),
              median_mean = median(changing_IF_mean))
}

plot_IF <- function(tpm_changing_IF){
  tpm_medians <- tpm_changing_IF %>% 
    summarise(median_IF_baseline = median(median_IF_baseline),
              mean_baseline = mean(mean_baseline),
              median_condition = median(median_condition),
              mean_condition = mean(mean_condition),
              median_mean = median(median_mean))
  
  plt_IF <- tpm_changing_IF %>% 
    ggplot() +
    aes(x = median_mean) +
    geom_density(alpha = 0.5, fill = "darkgray") +
    geom_vline(xintercept = tpm_medians$median_mean, linetype = "dashed") +
    labs(x = paste0("Fraction of \"Both\" genes' expression\nmediated by differential splicing",
                    "\n(Median: ", round(tpm_medians$median_mean, 3), ")"),
         y = "Density") +
    coord_cartesian(xlim = c(0, 1)) +
    theme(legend.position = "none")
  return(plt_IF)
}


if(fig2){
  
  
  plt_switches <- plot_switch_fraction(dgs_isoforms)
  # plt_switches
  ggsave(paste0("figs/isoform_switches.", fig_format), plt_switches)
  
  
  tpm_changing_IF <- compute_IF(tpms, dgs_isoforms, both_genes)
  
  plt_IF <- plot_IF(tpm_changing_IF)
  
  ggsave(paste0("figs/isoform_switches_IF.", fig_format), plt_switches/plt_IF)
  
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


if(fig2){
  
  switched_majorities <- count_switched_majority(dgs_isoforms, tpms)
  
  plt_switch_majority <- plot_switch_majority(switched_majorities)
  # plt_switch_majority
  ggsave(paste0("figs/switch_majority.", fig_format), plt_switch_majority)
}


## Full figure ----

if(fig2){
  
  layout <- "
AABB
AABB
CCDD
EEFF
"
  theme_set(theme_classic(base_size = 8))
  Figure2 <- plt_gene_counts + plt_gene_fractions + dgs_affected+plt_switches+plt_IF+dgs_signal +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') & 
    full_figure_theme()
  # Figure2
  
  ggsave(paste0("figs/Figure2.", fig_format), Figure2)
  
  theme_set(theme_classic(base_size = 19))
}


# Figure 3 ----

## Figure 3A: Significant pathways per experiment

plot_pathway_count <- function(ora_all){
  
  found_pathways <- ora_all %>% 
    group_by(experiment) %>% 
    summarise(DGE = sum(padj_deseq < 0.05, na.rm = TRUE),
              DGS = sum(padj_dexseq < 0.05, na.rm = TRUE),
              Overlap = sum((padj_dexseq < 0.05) & (padj_deseq < 0.05), na.rm = TRUE),
              Only_DGS = sum(!(padj_deseq < 0.05) & (padj_dexseq < 0.05), na.rm = TRUE)
    )
  
 # Sina plot of found pathways
  plt_pathway_counts <- found_pathways %>% 
    pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Pathways") %>% 
    mutate(Analysis = factor(Analysis, levels = c("DGE",
                                                  "DGS",
                                                  "Overlap",
                                                  "Only_DGS")
    ),
    Analysis =  recode(Analysis,
                       DGE = paste0("Differential\nExpression\n", "(Median: ", median(found_pathways$DGE), ")"),
                       DGS = paste0("Differential\nSplicing\n", "(Median: ", median(found_pathways$DGS), ")"),
                       Overlap = paste0("Overlap\n", "(Median: ", median(found_pathways$Overlap), ")"),
                       Only_DGS = paste0("Only\nDifferential\nSplicing\n", "(Median: ", median(found_pathways$Only_DGS), ")"))
    ) %>% 
    ggplot() +
    aes(y = Pathways, x = Analysis, color = Analysis) +
    ggforce::geom_sina() +
    scale_y_log10() +
    geom_boxplot(width = 0.05) +
    color_scale(reorder = c(2, 3, 1, 4)) + 
    labs(x = "",
         y = "Significant pathways\nper comparison") +
    theme(legend.position = "none")
  
  return(plt_pathway_counts)
  
}

if(fig3){
  
  ora_all <- readRDS("results/ora_all.RDS")
  
  plt_pathway_count <- plot_pathway_count(ora_all)
  # plt_pathway_count
  ggsave(plot = plt_pathway_count, filename = paste("figs/pathway_count.", fig_format))
}

## Figure 3B: Telomere pathways ----

if(fig3){
  # ora_all <- readRDS("results/ora_all.RDS")
  exp <- "77_GSE61220_TNF Treatment 12hrs"
  ora <- ora_all %>% 
    filter(experiment == exp)
  
  
  telomere_pathways <- plot_ora(ora, pattern = "Telomer", colors = c("darkgray", "purple", "lightblue"))
  ggsave(plot = telomere_pathways, filename = paste("figs/telomere_pathways.", fig_format))
}

## Figure 3C: Repair pathways----

if(fig3){
  # ora_all <- readRDS("results/ora_all.RDS")
  exp <- "79_GSE139262_SMARCB1 overexpression"
  ora <- ora_all %>% 
    filter(experiment == exp)
  
  
  repair_pathways <- plot_ora(ora, pattern = "Repair", colors = c("darkgray", "purple", "lightblue"), plotly=T)
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
    labs(x = paste0("Spearman's \u03C1 between gene-set \n enrichment scores", "\n(Median: ", round(median_correlation, 3), ")"),
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
    labs(x = "Spearman's \u03C1 between gene-set\nenrichment scores",
         y = "Count") +
    facet_grid(~association, space="free_x")  +
    fill_scale() +
    color_scale() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle=-90))

  return(plt_correlation)
}

if(fig3){
  
  # ora_all <- readRDS("results/ora_all.RDS")
  
  ora_correlation <- plot_ora_correlation(ora_all)
  # ora_correlation
  ggsave(plot = ora_correlation, filename = paste0("figs/ora_correlation.", fig_format))
  ora_correlation_facet <- plot_ora_correlation_facet(ora_all, threshold = 0)
  # ora_correlation_facet
  ggsave(plot = ora_correlation_facet, filename = paste0("figs/ora_correlation_facet.", fig_format), width = 6, height = 4)
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


if(fig3){
  
  # ora_all <- readRDS("results/ora_all.RDS")
  rr_shifts <- compute_rr_shifts(ora_all)
  
  
  rr_ridges <- plot_rr_ridges(rr_shifts)
  # rr_ridges
  ggsave(plot = rr_ridges, filename = paste0("figs/rr_ridges.", fig_format))
  rr_median <- plot_rr_median(rr_shifts)
  # rr_median
  ggsave(plot = rr_median, filename = paste0("figs/rr_median.", fig_format))
}



## Full figure ----

if(fig3){
  
  # Overview
  theme_set(theme_classic(base_size = 8))
  
  
  # layout <- "
  # #AAAA#
  # BBBCCC
  # DDDEEE
  # "
  layout <- "
  #AAAA#
  #BBBB#
  DDDEEE
  "
 
  Figure3 <- plt_pathway_count + telomere_pathways + #repair_pathways + 
    ora_correlation + rr_median +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') &
    full_figure_theme()
  # Figure3
  
  ggsave(paste0("figs/Figure3.", fig_format), Figure3, width = 8, height = 7)
  
  theme_set(theme_classic(base_size = 19))
}
