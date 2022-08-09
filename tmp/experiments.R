### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()
library(tidyverse)

source("tmp/run_experiment.R")
source("tmp/run_analysis.R")
source("tmp/concat.R")

rm(aggregate_pvalue)
theme_set(theme_classic(base_size = 20))
### Combine each and extract the comparisons to be run
experiments <- combine_experiments()

### Load GTF file
gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# row <- experiments[19,]

### Run experiments ----

### Define file to read from and group column
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
# groupCol <- "group_nr"

# apply(experiments[82,], 1, runExperiment, archs4db)
# row <- experiments[19,]
row <- experiments[1,]
### Analyse experiment results ----

# Load MSigDB
gene_sets <- pairedGSEA::prepare_msigdb()

# apply(row, 1, run_experiment, archs4db)
# apply(row, 1, analyse_experiment)

apply(experiments[100:199, ], 1, run_experiment, archs4db)
apply(experiments[100:199, ], 1, run_analysis, gene_sets)
apply(experiments, 1, getDDS, archs4db)

#4_GSE156101_TP53 and TNKS1-2T Knockout  

### Check missing
test <- stringr::str_c(experiments$study, experiments$`comparison_title (empty_if_not_okay)`)

test[which(table(test)>1)]


named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

rerun <- c(59,60, 73, 74, 184, 185)
apply(experiments[rerun, ], 1, runExperiment)
apply(experiments[rerun, ], 1, analyseExperiment)

apply(experiments, 1, missingExperiment)




### Check fgsea results
row <- experiments[8, ]
### Load metadata
md_file <- row$filename
dataname <- basename(md_file) %>% 
  stringr::str_remove(".xlsx") %>% 
  stringr::str_remove("csv") 
### Define experiment details
comparison <- row$`comparison (baseline_v_condition)`
experimentTitle <- row$`comparison_title (empty_if_not_okay)`
message("Analysing ", row$study, " ", experimentTitle)
### Results
comb <- readRDS(paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
fgseaRes <- readRDS(paste0("results/", dataname, "_fgseaRes_", experimentTitle, ".RDS"))
fgseaDxr <- readRDS(paste0("results/", dataname, "_fgseaDxr_", experimentTitle, ".RDS"))
fgseaDxr2 <- readRDS(paste0("results/", dataname, "_fgseaDxr2_", experimentTitle, ".RDS"))

### Significant gene sets
pathRes <- fgseaRes %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
pathDxr <- fgseaDxr %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
pathDxr2 <- fgseaDxr2 %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble

intersect(pathRes, pathDxr)
intersect(pathRes, pathDxr2)
intersect(pathDxr, pathDxr2)
intersect(pathRes, pathDxr) %>% intersect(pathDxr2)
union(pathRes, pathDxr) %>% union(pathDxr2)
union(pathDxr, pathDxr2) %>% intersect(pathRes)
setdiff(pathDxr, pathDxr2)
intersect(pathDxr, pathDxr2)

(overlap <- union(pathRes, pathDxr) %>% union(pathDxr2) %>% 
  mutate(deseq2 = 0 + (pathway %in% pathRes$pathway),
         dexseq = 0 + (pathway %in% pathDxr$pathway),
         dexseq2 = 0 + (pathway %in% pathDxr2$pathway),
         overlap = deseq2 + dexseq + dexseq2
         ) %>% 
  arrange(desc(overlap)))

library(ggplot2)
overlap %>% 
  ggplot(aes(x = pathway, y = 1, fill = overlap)) +
  geom_tile()

# concatFgsea <- concatFgseaResults(experiments)
concatFgsea <- readRDS("results/concatFgsea.RDS")

p <- concatFgsea %>% 
  mutate(overlap = as.factor(overlap)) %>% 
  ggplot(aes(x = pathway, y = experiment, fill = overlap)) +
  geom_tile() +
  scale_fill_manual(values = c("gray", "blue", "darkred"), na.value = "white")

ggsave("results/fgseaPlot.pdf", plot =  p, width = 200, height = 100, units = "cm", limitsize = FALSE)
p <- concatFgsea %>% 
  filter(overlap > 1) %>% 
  mutate(overlap = as.factor(overlap)) %>% 
  ggplot(aes(x = pathway, y = experiment, fill = overlap)) +
  geom_tile() +
  scale_fill_manual(values = c("blue", "darkred"), na.value = "white")
ggsave("results/fgseaPlot2.pdf", plot =  p, width = 200, height = 100, units = "cm", limitsize = FALSE)
p <- concatFgsea %>% 
  filter(dexseq2 == 1 & deseq2 == 0 & dexseq == 0) %>%
  mutate(dexseq2 = as.factor(dexseq2)) %>% 
  ggplot(aes(x = pathway, y = experiment, fill = dexseq2)) +
  geom_tile() 
  # scale_fill_manual(values = c("white", "darkblue"), na.value = "white")
ggsave("results/dexseq2unique.pdf", plot =  p, width = 200, height = 100, units = "cm", limitsize = FALSE)
p <- concatFgsea %>% 
  filter(deseq2 == 0 & dexseq == 1) %>%
  mutate(dexseq2 = as.factor(dexseq2)) %>% 
  ggplot(aes(x = pathway, y = experiment, fill = dexseq2)) +
  geom_tile() 
# scale_fill_manual(values = c("white", "darkblue"), na.value = "white")
ggsave("results/dexsequnique.pdf", plot =  p, width = 200, height = 100, units = "cm", limitsize = FALSE)


# Heatmap of representation
library(tidyverse)
p <- concatFgsea %>% 
  # filter(experiment == experiment[[1]]) %>% 
  pivot_longer(cols = c("deseq2", "dexseq", "dexseq2"),names_to = "analysis", values_to = "present") %>% 
  mutate(present = as.factor(present)) %>% 
  ggplot(aes(x = analysis, y = pathway, fill = present)) +
  geom_tile()
ggsave("results/present.pdf", plot =  p, width = 100, height = 300, units = "cm", limitsize = FALSE)

jac <- concatFgsea %>%
  mutate(across(starts_with("de"), .fns = ~as.logical(.x))) %>% 
  group_by(experiment) %>% 
  group_map(.f = function(df, ...){
    d <- df %>% 
      select(starts_with("de")) %>% 
      proxy::simil(by_rows = FALSE, method = "Jaccard")
    tibble(dedex = d[1], dedex2 = d[2], dexdex2 = d[3])
  }) %>% 
  bind_rows() %>% 
  mutate(experiment = concatFgsea %>% distinct(experiment) %>% pull)


sim <- concatFgsea %>%
  mutate(across(starts_with("de"), .fns = ~as.logical(.x))) %>% 
  group_by(experiment) %>% 
  group_map(.f = function(df, ...){
    d <- df %>% 
      select(starts_with("de")) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    tibble(dedex = d[1], dedex2 = d[2], dexdex2 = d[3])
  }) %>% 
  bind_rows() %>% 
  mutate(experiment = concatFgsea %>% distinct(experiment) %>% pull)

qplot(jac$dedex, main = "DESeq2 & DEXSeq overlap")
qplot(jac$dedex2, main = "DESeq2 & DEXSeq2 overlap")
qplot(jac$dexdex2, main = "DEXSeq & DEXSeq2 overlap")


## Study investigation

# Find relevant experiments
(relexp <- concatFgsea %>% 
  mutate(across(starts_with("de"), .fns = ~as.logical(.x))) %>% 
  filter(!deseq2 & (dexseq | dexseq2)) %>% 
  count(experiment) %>% 
  arrange(desc(n)))
(study <- concatFgsea %>%
  filter(experiment == relexp$experiment[1]))
study %>% 
  pivot_longer(cols = c("deseq2", "dexseq", "dexseq2"),names_to = "analysis", values_to = "present") %>% 
  mutate(present = as.factor(present)) %>% 
  ggplot(aes(x = analysis, y = pathway, fill = present)) +
  geom_tile()
median((select(concatFgsea, starts_with("de"))))

concatFgsea %>% 
  group_by(experiment) %>% 
  summarise(n1 = sum(deseq2), n2 = sum(dexseq), n3 = sum(dexseq2)) %>% 
  summarise(med1 = median(n1), med2 = median(n2), med3 = median(n3))



### ora analysis ----

ora_all <- concatenate_ora(experiments)
concatora <- readRDS("results/concatora.RDS")
ora_all <- readRDS("results/ora_all.RDS")
ora_all <- ora_all %>% 
  dplyr::select(-dplyr::starts_with("overlapG"))
concatora %>% 
  # filter(experiment == experiment[[1]]) %>% 
  pivot_longer(cols = starts_with("de"), names_to = "analysis", values_to = "present") %>% 
  mutate(present = as.factor(present),
         id = 1:nrow(.)) %>% 
  ggplot(aes(x = analysis, y = id, fill = present)) +
  geom_tile()


concatora %>% 
  filter(decombined & !dexseq & !deseq2) %>% 
  count(experiment)


Study <- "GSE101968 Obstructed defecation symdrome"
concatora %>% 
  filter(experiment == Study,
         decombined & !dexseq & !deseq2)

concatora %>% 
  filter(!decombined & dexseq & !deseq2) %>% 
  count(experiment)
concatora %>% 
  filter(!decombined & !dexseq & deseq2) %>% 
  count(experiment)

# ora Similarities ----
concatora %>% 
  group_by(experiment) %>% 
  summarise(jacdedex = proxy::simil(deseq2, dexseq, by_rows = FALSE, pairwise = TRUE, method = "Jaccard"),
            jacdecom = proxy::simil(deseq2, decombined, by_rows = FALSE, pairwise = TRUE, method = "Jaccard"),
            jacdexcom = proxy::simil(dexseq, decombined, by_rows = FALSE, pairwise = TRUE, method = "Jaccard"),
            simdedex = proxy::simil(deseq2, dexseq, by_rows = FALSE, pairwise = TRUE, method = "Simpson"),
            simdecom = proxy::simil(deseq2, decombined, by_rows = FALSE, pairwise = TRUE, method = "Simpson"),
            simdexcom = proxy::simil(dexseq, decombined, by_rows = FALSE, pairwise = TRUE, method = "Simpson")
            )

# Pathways per analysis
concatora %>% 
  group_by(experiment) %>% 
  summarise(deseq2 = sum(deseq2), dexseq = sum(dexseq), decombined = sum(decombined))

# Median pathways per analysis ----
ora_all %>% 
  group_by(experiment_title) %>% 
  summarise(deseq = sum(padj_deseq < 0.05), dexseq = sum(padj_dexseq < 0.05)) %>% 
  summarise(med_deseq2 = median(deseq), med_dexseq = median(dexseq))




# Median combined not in other
concatora %>% 
  filter((decombined & !dexseq & !deseq2) | (!decombined & (dexseq | dexseq))) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(decombined)) %>% 
  summarise(median = median(n))

# Median DEXSeq not in other
concatora %>% 
  filter((dexseq & !deseq2) | (!dexseq & !dexseq)) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(dexseq)) %>% 
  summarise(median = median(n))

# Median deseq2 not in other
concatora %>% 
  filter((!dexseq & deseq2) | (!dexseq & !dexseq)) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(deseq2)) %>% 
  summarise(median = median(n))

# Antallet af gene-sets hvor p-værdien er mindre end for hver af de individuelle analyser
oratot %>% 
  filter(across(starts_with("padj"), ~ .x < 0.05)) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(padj_decombined < padj_dexseq & padj_decombined < padj_deseq2))





## Rank shifts ----

# Plot distribution of adjusted p-values
oratot %>% 
  dplyr::filter(padj_deseq2 < 0.05 & padj_decombined < 0.05) %>% 
  tidyr::pivot_longer(cols = c("padj_deseq2", "padj_decombined"), names_to = "analysis", values_to = "padj") %>% 
  ggplot(aes(x = padj, fill = analysis))+
  geom_density(alpha = 0.5)

# Plot distribution of coverage (overlap / size)
oratot %>% 
  dplyr::filter(padj_deseq2 < 0.05 | padj_decombined < 0.05) %>% 
  tidyr::pivot_longer(cols = c("overlap_deseq2", "overlap_decombined"), names_to = "analysis", values_to = "overlap") %>% 
  dplyr::mutate(coverage = overlap / size) %>% 
  ggplot(aes(x = coverage, fill = analysis))+
  geom_density(alpha = 0.5)


# > oratot %>% filter(padj_decombined < 0.05) %>% nrow()
# [1] 377792
# > oratot %>% filter(padj_deseq2 < 0.05) %>% nrow()
# [1] 370800
# > oratot %>% filter(padj_decombined < 0.1) %>% nrow()
# [1] 462871
# > oratot %>% filter(padj_deseq2 < 0.1) %>% nrow()
# [1] 451905
# > oratot %>% filter(padj_deseq2 == 1) %>% nrow()
# [1] 120395
# > oratot %>% filter(padj_decombined == 1) %>% nrow()
# [1] 71002

# Rank each geneset within the experiments
shifts <- ora_all %>% 
  dplyr::group_by(experiment_title) %>% 
  dplyr::mutate(rank_deseq = rank(pval_deseq),
                rank_dexseq = rank(pval_dexseq),
                # rank_paired = rank(pval_paired),
                rank_shift = rank_deseq - rank_dexseq,
                # rank_shift_paired = rank_deseq - rank_paired,
                ES_shift = enrichment_score_deseq - enrichment_score_dexseq,
                # ES_shift_paired = enrichment_score_deseq - enrichment_score_paired)
                ) %>% 
  dplyr::filter(padj_deseq < 0.05 | padj_dexseq < 0.05) %>% 
  dplyr::arrange(desc(rank_shift))
# qplot(ranks$rank_shift)

# Store a list of all the top 20 ranked pathways of the paired analysis that were not in the top 20 for the deseq2 analysis
shifts %>% 
  dplyr::filter(
    rank_deseq > 50 & rank_paired < 50
  ) %>%  dplyr::select(pathway, experiment, dplyr::starts_with("rank"), dplyr::starts_with("odds_ratio"), dplyr::starts_with("enrichment_score"), dplyr::starts_with("ES_")) %>% #dplyr::arrange(experiment) %>% 
  dplyr::arrange(dplyr::desc(rank_shift_paired)) %>% 
  write.csv("results/ora_dge_v_paired_top50.csv")

# List of all the top 20 ranked pathways of the deseq2 analysis that were not in the top 20 for the paired analysis
ranks %>% 
  dplyr::filter(
    rank_deseq2 < 20 & rank_combined > 20
  ) %>%  dplyr::select(pathway, experiment, dplyr::starts_with("rank")) %>% dplyr::arrange(experiment)

shifts %>% 
  dplyr::filter(
    rank_deseq > 50 & rank_dexseq < 50
  ) %>%  dplyr::select(pathway, experiment_title, dplyr::starts_with("rank"), dplyr::starts_with("relative_risk"), dplyr::starts_with("enrichment_score"), ES_shift) %>% 
  dplyr::arrange(dplyr::desc(rank_shift)) %>% 
  write.csv("results/ora_dge_v_dtu_top50.csv")

readr::read_csv("~/Downloads/rstudio-export-2/ora_dge_v_dtu_top50.csv") %>% 
  select(-`...1`) %>% 
  writexl::write_xlsx("~/Downloads/rstudio-export-2/ora_dge_v_dtu_top50.xlsx")

readr::read_csv("~/Downloads/rstudio-export-2/ora_dge_v_paired_top50.csv") %>% 
  select(-`...1`) %>% 
  writexl::write_xlsx("~/Downloads/rstudio-export-2/ora_dge_v_paired_top50.xlsx")

readr::read_csv("~/Downloads/ora_dge_v_dtu_top50.csv") %>% 
  dplyr::select(-`...1`) %>% 
  writexl::write_xlsx("~/Downloads/ora_dge_v_dtu_top50.xlsx")

### Plotting ----

#' Hvad betyder det, at man kan finde flere go-termer med DTU og er de relevante?
#' Man skal OGSÅ kigger på DTU - ikke kun.
#' Et eller to grundige eksempel hvor vi ved hvad vi leder efter
#' Lav et par stykker til supplementary
#' Fig 1A: Præsentation af data - et diverst rum (UMAP)
#' Hierachical clustering viser hvad der er meget ens.
#' Gene-level analyse:
#'  (hvor mange er signifikante, overlap og forskelle) -
#'  Number DE/DU
#'  distribution af simpson coefficienter
#' Det samme på gensæt
#'  Number of genesets
#'  overlap
#'  Syngeri:
#'   Hvad er ændring i pværdier / onrichment score
#'   Slå op hvordan man beregner enrichment score / odds ratio
#'   Shift i pværdi / enrichment score
#' Cases:
#'  Træk gensæt skift i rank ud
#'  Hvilke gensæt nævner artiklen
#'  Find mening ud fra gensætnavne
#'  Kig primært på rank og rank shift
#'  Gensæt der kun er signifikante ved synergi
#'  
#' En reviewer kan måske spørge om man får noget ud af DTU alene - hav med i tabel  


concatResults <- concatRes(experiments)

concatResults <- readRDS("results/concatResults.RDS")

library(magrittr)
umap_data2 <- concatResults %>% 
  dplyr::select(transcript, experiment, log2FC) %>% 
  tidyr::pivot_wider(names_from = "experiment", values_from = "log2FC")
umap_data <- concatResults %>% 
  dplyr::select(transcript, experiment, log2FC) %>% 
  tidyr::pivot_wider(names_from = "transcript", values_from = "log2FC")

um <- umap_data %>% 
  dplyr::select(-experiment) %>% 
  dplyr::filter(across(everything(), ~ !is.na(.x)))
  uwot::umap(min_dist = 0.2)

um <- umap_data2 %>% 
    dplyr::select(-transcript) %>% 
    dplyr::filter(across(everything(), ~ !is.na(.x))) %>% 
    as.matrix() %>% 
    t() %>% 
    uwot::umap(n_neighbors = 3)

saveRDS(um, "results/all_umap2.RDS")

p <- um %>% 
  as_tibble(rownames = "Experiment") %>% 
  ggplot(aes(x = V1, y = V2))+
  geom_point() + 
  labs(title = "n_neighbors = 3")
  
# Zero significant genes
concatResults %>% 
  group_by(experiment) %>% 
  mutate(sig = padj_deseq2 < 0.05) %>% 
  summarise(n = sum(sig)) %>% 
  filter(n == 0)
"GSE147851 Cell-line difference in SMARCB1 overexpression background"







#  genes per dataset ----

concatenated_genes <- readRDS("results/concatenated_genes.RDS")

# DGS Genes per total ----

tot_genes <- concatenated_genes %>% filter(!is.na(pvalue_dexseq)) %>% distinct(gene) %>% nrow()

dgs_genes <- concatenated_genes %>% 
  filter(padj_dexseq < 0.05) %>% 
  mutate(dataset = str_extract(experiment, "\\d+_GSE\\d+")) %>% 
  distinct(gene, dataset, .keep_all = TRUE) %>% 
  dplyr::count(gene)

nrow(dgs_genes) / tot_genes

(dgs_genes %>% filter(n > 1) %>% nrow()) / tot_genes

# Median genes and recurrent genes per dataset

recurrent_genes <- concatenated_genes %>% 
  filter(padj_dexseq < 0.05) %>% 
  mutate(dataset = str_extract(experiment, "\\d+_GSE\\d+")) %>% 
  dplyr::count(dataset, gene) %>% 
  group_by(dataset) %>% 
  summarise(dexseq_recurrent_genes = sum(n > 1))

dexseq_genes <- concatenated_genes %>% 
  mutate(dataset = str_extract(experiment, "\\d+_GSE\\d+")) %>% 
  group_by(dataset) %>% 
  summarise(dexseq_genes = sum(padj_dexseq < 0.05, na.rm = TRUE),
            num_genes = length(unique(gene))) %>% 
  left_join(recurrent_genes, by = "dataset")

dexseq_genes %>% 
  filter(dexseq_recurrent_genes > 0) %>% 
  summarise(med_dexseq_recurrent = median(dexseq_recurrent_genes),
            frac_dexeq_recurrent = median(dexseq_recurrent_genes/num_genes)
  )


dexseq_genes %>% 
  pivot_longer(starts_with("dexseq"), names_to = "var", values_to = "gene_count") %>% 
  ggplot(aes(x = gene_count, fill = var)) +
  geom_density(alpha = 0.5) +
  # geom_vline(xintercept = median(gene_similarity$Similarity, na.rm = TRUE), linetype = "dashed") +
  # geom_histogram(binwidth = 0.02) +
  # labs(#title = "Overlap coefficient between DGE and DGS significant genes",
  #   x = "Overlap coefficient",
  #   y = "Density") +
  # scale_fill_manual(values = c("gray") + 
  # coord_cartesian(xlim = c(0,1)) +
  theme_classic(base_size = 20) 
  # theme(legend.position = "none")
  
# SVA metrics -------------------------------------------------------------

# concatenated_design <- concatenate_design(experiments)
# concatenated_design <- readRDS("results/concatenated_design.RDS")

(total_svas <- concatenated_design %>% 
  dplyr::mutate(design = as.character(design),
                svs = stringr::str_count(design, "sv\\d")) %>% 
  dplyr::count(svs))




# Gene OC ------------------------------------------------------------------


fig1c <- gene_similarity %>% 
  ggplot(aes(x = Similarity, fill = "Gray")) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = median(gene_similarity$Similarity, na.rm = TRUE), linetype = "dashed") +
  # geom_histogram(binwidth = 0.02) +
  labs(#title = "Overlap coefficient between DGE and DGS significant genes",
       x = "Overlap coefficient",
       y = "Density") +
  scale_fill_manual(values = "gray") + 
  coord_cartesian(xlim = c(0,1)) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

ggsave(plot = fig1c, filename = "figs/fig1c.png")

# Transcript fractions ------------------------------------------------------------------

# Distribution of fraction of tested transcripts that are DTE? 
# (with the notion if it is only 1 of 10 transcripts that are DTE 
# it is not the gene that is differentially expressed but a single isoofrm.
# Vise verca if all isoforms are upregulated it is DGE not DGU).

concatenated_genes <- readRDS("results/concatenated_genes.RDS")

significant_genes <- concatenated_genes %>% 
  filter(padj_dexseq < 0.05 | padj_deseq < 0.05)
  

# concatenated_results <- concatenate_results(experiments)
concatenated_results <- readRDS("results/concatenated_results.RDS")

# analysed_genes <- concatenated_results %>% 
#   count(experiment, gene) %>% 
#   filter(n > 1) %>% 
#   distinct(experiment, gene)

transcript_fractions <- concatenated_results %>% 
  semi_join(significant_genes, by = c("experiment", "gene")) %>%
  group_by(experiment, gene) %>% 
  summarise(fraction_dtu = sum(padj_dexseq < 0.05, na.rm = TRUE) / n(), 
            fraction_dte = sum(padj_deseq < 0.05, na.rm = TRUE) / n(),
            n = n()) %>% 
  filter(n > 1) %>%
  left_join(significant_genes, by = c("experiment", "gene")) %>% 
  mutate(association = dplyr::case_when(
    padj_dexseq < 0.05 & padj_deseq < 0.05 ~ "Both",
    padj_dexseq < 0.05 ~ "DGS",
    padj_deseq < 0.05 ~ "DGE",
    TRUE ~ "NA"
  ))

saveRDS(transcript_fractions, "results/transcript_fractions.RDS")

transcript_fractions <- readRDS("results/transcript_fractions.RDS")

exclude <- transcript_fractions %>% 
  count(experiment, association) %>% 
  filter(n < 10)

# densities <- transcript_fractions %>% 
#   anti_join(exclude, by = "experiment") %>%
#   group_by(experiment) %>%
#   summarise(dens = graphics::hist(fraction, breaks = seq(-0.06, 1.06, by = 0.12),
#                plot = FALSE)$density,
#             frac = graphics::hist(fraction, breaks = seq(-0.06, 1.06, by = 0.12),
#                                   plot = FALSE)$mids)

# This one!!
densities <- transcript_fractions %>% 
  anti_join(exclude, by = c("experiment", "association")) %>%
  group_by(experiment, association) %>%
  summarise(dens = density(fraction_dte, from = 0, to = 1)$y,
            frac = density(fraction_dte, from = 0, to = 1)$x,
            association = association[1])

fig1d <- densities %>% 
  filter(association != "NA") %>% 
  ggplot() +
  aes(x = frac, y = dens, color = experiment) +
  # geom_smooth(method = "loess", se = FALSE) +
  geom_line() +
  geom_smooth(color = "black", se = FALSE, size = 1) +
  scale_color_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
  labs(y = "Density",
       x = "Fraction of significant transcripts in a gene") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  facet_wrap(~association)

ggsave(plot = fig1d, filename = "figs/fig1d.png")

# Median fraction per study

transcript_fractions %>% 
  anti_join(exclude, by = c("experiment", "association")) %>%
  filter(association != "DGS") %>%
  mutate(n = case_when(n == 2 ~ "2",
                       n == 3 ~ "3",
                       n %in% 4:6~ "4-6",
                       n >6 ~ "7+"),
         association = paste0(n, association)) %>% 
  group_by(experiment, association) %>% 
  summarise(mean = mean(fraction_dte), med = median(fraction_dte), med_dtu = median(fraction_dtu)) %>% 
  # box_data %>% 
  # filter()
  ggplot() +
  aes(x = association, y = med) +
  geom_boxplot() +
  theme_classic(base_size = 20)

  
# Densities split
densities <- transcript_fractions %>% 
  anti_join(exclude, by = c("experiment", "association")) %>%
  # filter(association != "DGS") %>% 
  mutate(n = case_when(n == 2 ~ "2",
                       n == 3 ~ "3",
                       n %in% 4:6~ "4-6",
                       n > 6 ~ "7+"),
         fraction = case_when(association == "DGS" ~  fraction_dtu,
                              TRUE ~ fraction_dte)
         # association = paste0(n, association)
         ) %>% 
  group_by(experiment, association, n) %>%
  filter(n() > 2) %>%
  summarise(dens = density(fraction_dtu, from = 0, to = 1)$y,
            frac = density(fraction_dtu, from = 0, to = 1)$x,
            association = association[1],
            n = n[1])

densities %>% 
  filter(association != "NA") %>% 
  ggplot() +
  aes(x = frac, y = dens, color = association) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_line() +
  geom_smooth(se = FALSE, size = 1) + #color = "black"
  # scale_color_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
  labs(y = "Density",
       x = "Fraction of significant transcripts in a gene") +
  theme_classic(base_size = 20) +
  # theme(legend.position = "none") +
  facet_wrap(~n)


transcript_fractions %>% 
  anti_join(exclude, by = c("experiment", "association")) %>%
  filter(!is.na(association) & !is.na(n)) %>%
  mutate(n = case_when(n == 2 ~ "2",
                       n == 3 ~ "3",
                       n %in% 4:6~ "4-6",
                       n > 6 ~ "7+"),
         fraction = case_when(association == "DGS" ~  fraction_dtu,
                              TRUE ~ fraction_dte)
         # association = paste0(n, association)
         ) %>% 
  group_by(experiment, association, n) %>%
  filter(n() > 2) %>%
  summarise(mean = mean(fraction_dtu), med = median(fraction_dtu), n = n) %>% 
  # box_data %>% 
  # filter()
  ggplot() +
  aes(x = association, y = med) +
  geom_boxplot() +
  facet_wrap(~n) +
  theme_classic(base_size = 20)


## TF V2 ----
if(TRUE){

n_transcripts <- transcript_fractions %>% 
  # anti_join(exclude, by = c("experiment", "association")) %>%
  filter(!is.na(association) & !is.na(n)) %>%
  mutate(n_c = case_when(n > 3 ~ "4+",
                         TRUE ~ as.character(n)),
         transcripts = case_when(fraction_dte*n <= 3 ~ as.character(fraction_dte*n),
                                 TRUE ~ "4+"),
         transcripts_dtu = case_when(fraction_dtu*n <= 3 ~ as.character(fraction_dtu*n),
                                     TRUE ~ "4+")
         # association = paste0(n, association)
  )

gene_counts <- dplyr::count(n_transcripts, experiment, association, n_c, name = "gene_count")
n_counts <- dplyr::count(n_transcripts, experiment, association, n_c, transcripts, name = "n_count")
n_counts_dtu <- dplyr::count(n_transcripts, experiment, association, n_c, transcripts_dtu, name = "n_count")

# Boxplot
# n_transcripts %>% 
#   left_join(gene_counts, by = c("experiment", "association", "n_c")) %>% 
#   left_join(n_counts,  by = c("experiment", "association", "n_c", "transcripts")) %>% 
#   mutate(gene_fraction = n_count/gene_count) %>% 
#   ggplot() +
#   aes(x = transcripts, color = association, y = gene_fraction) +
#   geom_boxplot() +
#   facet_wrap(~n_c) +
#   theme_classic(base_size = 20)

# DGE
fig2c_dge <- n_transcripts %>% 
  left_join(gene_counts, by = c("experiment", "association", "n_c")) %>% 
  left_join(n_counts,  by = c("experiment", "association", "n_c", "transcripts")) %>% 
  mutate(gene_fraction = n_count/gene_count) %>% 
  group_by(association, n_c, transcripts) %>% 
  summarise(med_gene_fraction = median(gene_fraction),
            association = association[1], n_c = n_c[1], transcripts = transcripts[1]) %>% 
  ggplot() +
  aes(x = transcripts, color = association, y = med_gene_fraction, group = association) +
  # geom_boxplot() +
  geom_line() +
  geom_point(size = 3) +
  # facet_wrap(~n_c) +
  facet_grid(~n_c, space="free_x")  +
  labs(x = "Differentially expressed transcripts in gene",
       y = "Median fraction of genes",
       color = "padj < 0.05") +
  scale_color_manual(values = c("gray", "lightblue", "darkred")) +
  theme_classic(base_size = 20)

ggsave(plot = fig2c_dge, filename = "figs/fig2c_dge.png")



fig2c_dgs <- n_transcripts %>% 
  left_join(gene_counts, by = c("experiment", "association", "n_c")) %>% 
  left_join(n_counts_dtu,  by = c("experiment", "association", "n_c", "transcripts_dtu")) %>% 
  filter(n_c == "3") %>% 
  mutate(gene_fraction = n_count/gene_count) %>% 
  group_by(association, n_c, transcripts_dtu) %>% 
  summarise(med_gene_fraction = median(gene_fraction),
            association = association[1], n_c = n_c[1], transcripts = transcripts_dtu[1]) %>% 
  ggplot() +
  aes(x = transcripts, color = association, y = med_gene_fraction, group = association) +
  geom_line() +
  geom_point(size = 3) +
  # facet_grid(~n_c, space="free_x")  +
  labs(x = "Differentially spliced transcripts in gene",
       y = "Median fraction of genes",
       color = "padj < 0.05") +
  scale_color_manual(values = c("gray", "lightblue", "darkred")) +
  theme_classic(base_size = 20)

ggsave(plot = fig2c_dgs, filename = "figs/fig2c_dgs.png")

if(FALSE){ # Boxplot versions
  # DGE
  n_transcripts %>% 
    left_join(gene_counts, by = c("experiment", "association", "n_c")) %>% 
    left_join(n_counts,  by = c("experiment", "association", "n_c", "transcripts")) %>% 
    mutate(gene_fraction = n_count/gene_count) %>% 
    ggplot() +
    aes(x = transcripts, color = association, y = gene_fraction) +
    geom_boxplot() +
    facet_grid(~n_c, space="free_x")  +
    labs(x = "Differentially expressed transcripts in gene",
         y = "Median fraction of genes",
         color = "padj < 0.05") +
    scale_color_manual(values = c("gray", "lightblue", "darkred")) +
    theme_classic(base_size = 20)
  
  n_transcripts %>% 
    left_join(gene_counts, by = c("experiment", "association", "n_c")) %>% 
    left_join(n_counts_dtu,  by = c("experiment", "association", "n_c", "transcripts_dtu")) %>%
    mutate(gene_fraction = n_count/gene_count) %>% 
    ggplot() +
    aes(x = transcripts_dtu, color = association, y = gene_fraction) +
    geom_boxplot() +
    facet_grid(~n_c, space="free_x")  +
    labs(x = "Differentially spliced transcripts in gene",
         y = "Median fraction of genes",
         color = "padj < 0.05") +
    scale_color_manual(values = c("gray", "lightblue", "darkred")) +
    theme_classic(base_size = 20)
  
  
}


}
  
transcript_fractions %>% 
  anti_join(exclude, by = c("experiment", "association")) %>%
  filter(!is.na(association) & !is.na(n)) %>%
  mutate(n_c = case_when(n %in% 4:6~ "4-6",
                         n > 6 ~ "7+",
                         TRUE ~ as.character(n)),
         transcripts = case_when(fraction_dte*n < 6 ~ as.character(fraction_dte*n),
                                 TRUE ~ "7+")
         # association = paste0(n, association)
  ) %>% 
  group_by(experiment, association) %>% 
  filter(n() > 1) %>%
  # summarise(mean = mean(fraction_dtu), med = median(fraction_dtu), n = n) %>% 
  # box_data %>% 
  # filter()
  ggplot() +
  aes(x = transcripts, color = association, y = fraction_dte) +
  geom_boxplot() +
  facet_wrap(~n_c) +
  theme_classic(base_size = 20)
  
# not this!
if(FALSE){
  transcript_fractions %>% 
    anti_join(exclude, by = "experiment") %>%
    ggplot() +
    aes(x = fraction, y = ..density.., color = experiment) +
    # aes(x = frac, y = dens color = experiment) +
    # geom_smooth(method = "loess", se = FALSE) +
    geom_line(stat = "density") +
    # geom_line(stat = "density", color = "black") +
    # geom_smooth(data = densities, mapping = aes(x = frac, y = dens), color = "black", se = FALSE) +
    scale_color_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
    labs(y = "Density",
         x = "Fraction of significant transcripts in a gene") +
    theme_classic(base_size = 20) +
    theme(legend.position = "none")
  
  
  
  fig1d <- transcript_fractions %>% 
    anti_join(exclude, by = "experiment") %>%
    ggplot() +
    aes(x = fraction, fill = experiment) +
    geom_density(alpha = 0.02, color = NA) +
    geom_density(fill = NA) +
    scale_fill_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
    labs(y = "Density",
         x = "Fraction of significant transcripts in a gene") +
    theme_classic(base_size = 20) +
    theme(legend.position = "none")
  
  ggsave(plot = fig1d, filename = "figs/fig1d.png")
  
  
  # This one!
  fig1d <- transcript_fractions %>% 
    anti_join(exclude, by = "experiment") %>%
    ggplot() +
    aes(x = fraction, color = experiment) +
    geom_density(alpha = 0.03) +
    geom_density(bw = 1, color = "black") +
    scale_color_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
    labs(y = "Density",
         x = "Fraction of significant transcripts in a gene") +
    theme_classic(base_size = 20) +
    theme(legend.position = "none")
  
  ggsave(plot = fig1d, filename = "figs/fig1d.png")
  
  fig1d <- transcript_fractions %>% 
    # anti_join(exclude, by = "experiment") %>% 
    filter(gene != "NA") %>% 
    ggplot() +
    aes(x = fraction, y = n, color = experiment) +
    geom_point(alpha = 0.1) +
    # geom_density(fill = NA) +
    scale_color_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
    labs(y = "Number of transcripts in gene",
         x = "Fraction of significant transcripts in a gene") +
    theme_classic(base_size = 20) +
    theme(legend.position = "none")
  
}


# Pathway OC ------------------------------------------------------------------


pathway_similarity <- ora_all %>%
  mutate(DGE = padj_deseq < 0.05,
         DGS = padj_dexseq < 0.05,
         # paired = padj_paired < 0.05,
         added = padj_deseq < 0.05 | padj_dexseq < 0.05) %>% 
  group_by(experiment) %>% 
  summarise(Similarity = proxy::simil(x = DGE, y = DGS, 
                                      by_rows = FALSE, method = "Simpson")[1]) %>% 
  filter(!is.na(Similarity))

# as density
fig2b <- pathway_similarity %>% 
  # pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  ggplot(aes(x = Similarity, fill = "Gray")) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = median(pathway_similarity$Similarity),
             linetype = "dashed") +
  # scale_color_brewer(palette = "BuPu") +
  labs(x = "Overlap coefficient",
       y = "Density") +
  scale_color_manual(values = "gray") + 
  scale_fill_manual(values = "gray") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

ggsave(plot = fig2b, filename = "figs/fig2b.png")

# as sina
pathway_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  filter(!is.na("Similarity")) %>% 
  ggplot(aes(y = Similarity, x = Method)) +
  ggforce::geom_sina() +
  theme_classic()

# Enrichment density ------------------------------------------------------------------

# ora_all <- concatenate_ora(experiments)
ora_all <- readRDS("results/ora_all.RDS")


fig2e <- ora_all %>% 
  filter(padj_dexseq < 0.05) %>%
  mutate(enrichment_score_s = relative_risk_dexseq - relative_risk_deseq) %>% 
  filter(abs(enrichment_score_s) > 1) %>% 
  mutate(enrichment_score_s = log2(abs(enrichment_score_s))* sign(enrichment_score_s)) %>%
  ggplot() +
  aes(x = enrichment_score_s, fill = experiment) +
  geom_density(alpha = 0.02, color = NA) +
  geom_density(fill = NA) +
  scale_fill_manual(values = rep("gray", length(unique(ora_all$experiment)))) +
  # scale_x_continuous(trans = "log2") +
  # coord_cartesian(xlim = c(-4,10)) +
  labs(y = "Density",
       x = "Enrichment Score Shift") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

ggsave(plot = fig2e, filename = "figs/fig2e.png")



## fgsa analysis ----

# concatFgsea <- concatFgseaResults(experiments)

concatFgsea %>% 
  group_by(experiment) %>% 
  summarise(n1 = sum(deseq), n2 = sum(dexseq), n3 = sum(dexseqlfc)) %>% 
  summarise(med_deseq = median(n1), med_dexseq = median(n2), med_dexseqlfc = median(n3))



## RR shifts ----

rr_shifts <- ora_all %>% 
  dplyr::group_by(experiment) %>% 
  dplyr::mutate(relative_risk_shift = relative_risk_deseq - relative_risk_dexseq) %>% 
  dplyr::filter(padj_dexseq < 0.05) %>% 
  dplyr::select(pathway, experiment, dplyr::starts_with("relative")) %>% 
  dplyr::arrange(desc(relative_risk_shift))


# V2
rr_shifts %>% 
  dplyr::ungroup() %>%
  # dplyr::filter(relative_risk_shift < -2 | relative_risk_shift > 2) %>%
  dplyr::mutate(Shift = -relative_risk_shift,
                Shift = log2(abs(Shift))*sign(Shift),
                experiment = paste0(stringr::str_extract(experiment, "\\d+_GSE\\d+_\\w"), stringr::str_sub(experiment, -1))) %>% 
  dplyr::select(Shift, pathway, experiment) %>% 
  # pivot_longer(cols = -c("experiment", "odds_ratio_paired", "pathway"), names_to = "Analysis", values_to = "Pathways") %>% 
  # mutate(Analysis = factor(Analysis, levels = c("DESeq2", "DEXSeq", "Shift")
                           #c("DESeq2", "Paired", "Overlap", "Only paired", "Only DESeq2")
  # )) %>% 
  ggplot(aes(y = pathway, x = experiment, fill = Shift)) +
  # stat_summary(fun = median, geom = "point", size = 3, color = "red", shape = 18) +
  # geom_violin() +
  # ggforce::geom_sina() +
  geom_tile() +
  # geom_boxplot(width = 0.05) +
  labs(x = "") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low = "darkred", high = "navy")

# V3

rr_shifts <- ora_all %>% 
  dplyr::group_by(experiment) %>% 
  dplyr::summarise(median_rr_shift = median(relative_risk_deseq - relative_risk_dexseq)) 


rr_shifts %>% 
  ggplot(aes(x = median_rr_shift)) +
  geom_histogram(fill = "darkgray", color = "white")


# V4

rr_shifts <- ora_all %>% 
  dplyr::mutate(Shift = relative_risk_dexseq - relative_risk_deseq,
                association = dplyr::case_when(
                  padj_dexseq < 0.05 & padj_deseq < 0.05 ~ "Both",
                  padj_dexseq < 0.05 ~ "DGS",
                  padj_deseq < 0.05 ~ "DGE",
                  TRUE ~ "Neither"),
                association = factor(association, levels = c("Both", "DGE", "DGS", "Neither"))) %>% 
  dplyr::filter(association != "Neither") %>%
  dplyr::select(pathway, experiment, Shift, association) %>% 
  dplyr::arrange(desc(Shift))

rr_shifts %>% 
  dplyr::mutate(Shift = log2(abs(Shift))*sign(Shift)) %>%
  ggplot(aes(x = Shift, fill = association)) +
  geom_histogram(color = "white") +
  facet_grid(~association, space="free_x")  +
  scale_fill_manual(values = c("gray", "lightblue", "darkred"), na.value = "white") +
  labs(x = "Enrichment score shift",
       y = "Count",
       # fill = "padj < 0.05"
       ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

# V5
rr_shifts %>% 
  dplyr::mutate(Shift = log2(abs(Shift))) %>%
  ggplot(aes(y = Shift, x = association, color = association)) +
  ggforce::geom_sina() +
  geom_boxplot(width = 0.05) +
  # facet_grid(~association, space="free_x")  +
  scale_color_manual(values = c("gray", "lightblue", "darkred", "black"), na.value = "white") +
  labs(x = "padj < 0.05",
       y = "Enrichment score shift"
       # y = "Relative Risk shift"
       # fill = "padj < 0.05"
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

# V6
rr_shifts %>%
  group_by(experiment, association) %>% 
  summarise(Shift = median(abs(Shift))) %>% 
  # filter(Shift != 0) %>%
  # dplyr::mutate(Shift = abs(log2(Shift))) %>%
  ggplot(aes(x = Shift, y = association, color = association, fill = association)) +
  ggridges::geom_density_ridges(scale = 4, alpha = 0.4, quantile_lines = TRUE, quantiles = 2) +
  # facet_grid(~association, space="free_x")  +
  scale_color_manual(values = c("gray", "lightblue", "darkred", "black"), na.value = "white") +
  scale_fill_manual(values = c("gray", "lightblue", "darkred", "black"), na.value = "white") +
  labs(x = "Relative risk shift",
       y = ""
       # y = "Relative Risk shift"
       # fill = "padj < 0.05"
  ) +
  scale_x_log10() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

rr_shift_median <- rr_shifts %>%
  group_by(experiment) %>% 
  summarise(Shift = median(abs(Shift)))
rr_shift_median %>% 
  ggplot(aes(x = Shift)) +
  geom_histogram(fill = "darkgray", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(rr_shift_median$Shift), linetype = "dashed") +
  labs(x = paste0("Relative risk shift (Median: ", round(median(rr_shift_median$Shift),2), ")")
  ) +
  scale_x_log10() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

rr_shifts %>% 
  group_by(experiment, association) %>% 
  summarise(Shift = median(abs(Shift)+1)) %>% 
  filter(!is.na(Shift)) %>%
  # dplyr::mutate(Shift = abs(log2(Shift)))  %>%
  group_by(association) %>% 
  summarise(median = median(Shift))


# V7 (in plot_paper)



  
# SVA contrasts----

# Test
apply(experiments[1, ], 1, run_experiment, archs4db, deseq_only = T, run_sva = F)



concatenated_no_sva_genes <- deseq_no_sva(experiments, archs4db, gtf)

concatenated_no_sva_genes <- readRDS("results/concatenated_no_sva_genes.RDS")

concatenated_genes <- readRDS("results/concatenated_genes.RDS")
concatenated_sva_genes <- concatenated_genes %>% filter(padj_deseq < 0.05)


gene_counts <- concatenated_sva_genes %>% 
  dplyr::count(experiment)
median_gene_counts <- median(gene_counts$n)
median_gene_counts

gene_counts_no_sva <- concatenated_no_sva_genes %>% 
  dplyr::count(experiment)
median_gene_counts_no_sva <- median(gene_counts_no_sva$n)
median_gene_counts_no_sva

gene_diff <- deseq_genes %>% 
  full_join(deseq_genes_no_sva, by = c("experiment", "gene")) %>% 
  group_by(experiment) %>% 
  summarise(count = sum(!is.na(padj_deseq)), count_no_sva = sum(!is.na(padj))) %>% 
  mutate(diff = count_no_sva - count)

gene_diff %>% 
  ggplot() +
  aes(y = diff, x = "SVA genes") +
  geom_boxplot() +
  theme_classic()

# Similarity

gene_change_sva <- concatenated_sva_genes %>% 
  mutate(gene_change = paste0(sign(lfc_deseq))) %>% 
  select(gene_change, experiment, gene)

gene_change_no_sva <- concatenated_no_sva_genes %>% 
  mutate(gene_change = paste0(sign(lfc))) %>% 
  select(gene_change, experiment, gene)



gene_change_sva %>% inner_join(gene_change_no_sva, by = c("experiment", "gene"), suffix = c("_sva", "_no_sva")) %>% 
  group_by(experiment) %>% 
  summarise(diff = sum(gene_change_sva != gene_change_no_sva, na.rm=TRUE)) %>% 
  summarise(median = median(diff))
gene_change_sva %>% anti_join(gene_change_no_sva, by = c("experiment", "gene")) %>% count(experiment) %>% summarise(median = median(n))
gene_change_no_sva %>% anti_join(gene_change_sva, by = c("experiment", "gene")) %>% count(experiment) %>% summarise(median = median(n))


# Median difference gene to pathway ----


sig_genes <- concatenated_genes %>% 
  group_by(experiment) %>% 
  summarise(deseq = sum(padj_deseq < 0.05, na.rm = TRUE),
            dexseq = sum(padj_dexseq < 0.05, na.rm = TRUE))

sig_genesets <- ora_all %>% 
  group_by(experiment) %>% 
  summarise(deseq = sum(padj_deseq < 0.05, na.rm = TRUE),
            dexseq = sum(padj_dexseq < 0.05, na.rm = TRUE))

median(sig_genesets$deseq)
median(sig_genesets$dexseq)

sig <- sig_genes %>% 
  left_join(sig_genesets, by = c("experiment"), suffix = c("_gene", "_geneset")) %>% 
  summarise(DGE = median(deseq_geneset / deseq_gene, na.rm = TRUE),
         DGS = median(dexseq_geneset / dexseq_gene, na.rm = TRUE))


# DGS vs DGE per gene ----

single_transcript_genes <- concatenated_genes %>% 
  filter(!is.na(padj_dexseq)) %>% 
  distinct(gene)

concatenated_genes %>% 
  semi_join(single_transcript_genes, by = "gene") %>% 
  group_by(gene) %>% 
  summarise(DGS = sum(padj_dexseq < 0.05, na.rm = TRUE),
            DGE = sum(padj_deseq < 0.05, na.rm = TRUE)) %>% 
  # filter(DGS != 0 & DGE != 0) %>% 
  ggplot() +
  aes(x = DGE, y = DGS, fill = cut(..count.., c(0,10,50,100, 200, 500, 1000,Inf)), bins=30) + 
  geom_hex() +
  geom_abline(slope = 1) +
  # labs(fill = "bins") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")


ora_all %>% 
  # semi_join(single_transcript_genes, by = "gene") %>% 
  group_by(pathway) %>% 
  summarise(DGS = sum(padj_dexseq < 0.05, na.rm = TRUE),
            DGE = sum(padj_deseq < 0.05, na.rm = TRUE)) %>% 
  ggplot() +
  aes(x = DGE, y = DGS, fill = cut(..count.., c(0,10,50,100, 200, 500, 1000,Inf)), bins=30) +
  geom_hex() +
  geom_abline(slope = 1) +
  labs(title = "Gene Set") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")


# Both genes ----
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
tpms <- readRDS("results/tpms.RDS")

concatenated_genes <- readRDS("results/concatenated_genes.RDS")
both_genes <- concatenated_genes %>% 
  filter(padj_dexseq < 0.05 & padj_deseq < 0.05,
         gene != "NA")

dgs_isoforms <- readRDS("results/concatenated_results.RDS") %>% 
  semi_join(both_genes, by = c("gene", "experiment"))

# 2

tpm_two <- tpms %>% 
  group_by(gene, experiment) %>% 
  # filter(max(tpm_baseline) != 0 & max(tpm_condition) != 0) %>% 
  summarise(major_baseline = transcript[which(tpm_baseline == max(tpm_baseline))], 
            major_condition = transcript[which(tpm_condition == max(tpm_condition))])

tpm_two_change <- tpm_two %>% 
  group_by(experiment) %>% 
  summarise(num_genes = n(), major_change = sum(major_baseline != major_condition)) %>% 
  mutate(change_percent = major_change/num_genes)

median(tpm_two_change$major_change)
median(tpm_two_change$change_percent)


# maxes
two_test <- tpms %>% 
  group_by(gene, experiment) %>% 
  # filter(max(tpm_baseline) != 0 & max(tpm_condition) != 0) %>% 
  summarise(major_baseline = length(which(tpm_baseline == max(tpm_baseline))), 
            major_condition = length(which(tpm_condition == max(tpm_condition))))
two_test <- tpms %>% 
  group_by(gene, experiment) %>% 
  # filter(max(tpm_baseline) != 0 & max(tpm_condition) != 0) %>% 
  summarise(major_baseline = max(tpm_baseline), 
            major_condition = max(tpm_condition))
double_major <- two_test %>% filter(str_detect(major_baseline, " ") | str_detect(major_condition, " "))

one_two %>% 
  filter(gene %in% double_major$gene & experiment %in% double_major$experiment) %>% 
  select(tpm_baseline, tpm_condition, gene, experiment, transcript) %>% 
  arrange(experiment, gene)



# 3


dgs_isoforms <- readRDS("results/concatenated_results.RDS") %>% 
  semi_join(both_genes, by = c("gene", "experiment"))
  

tpm_threshold <- function(tpms, threshold, plot = FALSE){
  
  tpm_X <- tpms %>% 
    right_join(.GlobalEnv$dgs_isoforms, by = c("gene", "transcript", "experiment")) %>% 
    group_by(gene, experiment) %>% 
    filter(tpm_baseline >= threshold*sum(tpm_baseline) | tpm_condition >= threshold*sum(tpm_condition) & padj_dexseq < 0.05) %>% 
    distinct(gene, experiment) %>% 
    ungroup() %>% 
    count(experiment, name = "events") %>% 
    left_join(.GlobalEnv$both_genes %>% count(experiment, name = "num_genes"), by = c("experiment")) %>% 
    mutate(percent = events/num_genes)
  
  if(plot){
    
    plt <- tpm_X %>% 
      ggplot() +
      aes(x = percent) +
      geom_histogram(fill = "darkgray", color = "white", alpha = 0.8) +
      geom_vline(xintercept = median(tpm_X$percent), linetype = "dashed") +
      theme_classic(base_size = 20)
    return(plt)
    
  }
  
  message("Median Events: ", median(tpm_X$events), "\n",
               "Median Percent: ", median(tpm_X$percent))
  
  return(tpm_X)
  
}


tpm_three_10 <- tpm_threshold(tpms, 0.1, T)
tpm_three_25 <- tpm_threshold(tpms, 0.25)
tpm_three_50 <- tpm_threshold(tpms, 0.5)


  
# 5

tpm_dIF <- tpms %>% 
  right_join(dgs_isoforms, by = c("transcript", "gene", "experiment")) %>% 
  group_by(gene, experiment) %>% 
  mutate(IF_baseline = tpm_baseline / sum(tpm_baseline),
         IF_condition = tpm_condition / sum(tpm_condition),
         dIF = IF_condition - IF_baseline) %>% 
  filter(padj_dexseq < 0.05)

gene_counts <- both_genes %>% count(experiment, name = "num_genes")

dIF_threshold <- function(tpm_dIF, threshold){
  
  dIF_X <- tpm_dIF %>% 
    filter(abs(dIF) >= threshold) %>% 
    distinct(gene, experiment) %>% 
    ungroup() %>% 
    count(experiment, name = "events") %>% 
    left_join(.GlobalEnv$gene_counts, by = c("experiment")) %>% 
    mutate(percent = events/num_genes)
  
  message("Median Events: ", median(dIF_X$events), "\n",
          "Median Percent: ", median(dIF_X$percent))
  
  return(dIF_X)
}

dIF_10 <- dIF_threshold(tpm_dIF, 0.1)
dIF_25 <- dIF_threshold(tpm_dIF, 0.25)
dIF_33 <- dIF_threshold(tpm_dIF, 0.33)
dIF_50 <- dIF_threshold(tpm_dIF, 0.5)


tpm_dIF %>% 
  ggplot() +
  aes(x = IF_baseline) +
  geom_histogram(fill = "darkgray", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(tpm_dIF$IF_baseline), linetype = "dashed") +
  theme_classic(base_size = 20)


# 6

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
  labs(y = paste0("Fraction of gene expression\nmediated by differential splicing",
                  "\n(Median: ", round(tpm_medians$median_mean, 3), ")"),
       x = "") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(legend.position = "none",
        axis.title.y = element_text(angle = 0))


experiment_title <- "76_GSE183984_Cetuximab treatment of primary colorectal cancer"
# experiment_title <- "79_GSE139262_SMARCB1 overexpression"

dexseqres <- readRDS(paste0("results/", experiment_title, "_dexseqres.RDS"))
deseqres <- readRDS(paste0("results/", experiment_title, "_deseq2res.RDS"))
ora_dexseq <- readRDS(paste0("results/", experiment_title, "_ora_dexseq.RDS"))
ora_deseq <- readRDS(paste0("results/", experiment_title, "_ora_deseq.RDS"))
aggregated_pvals <- readRDS(paste0("results/", experiment_title, "_aggregated_pvals.RDS"))


## SVA Issue ----

dds <- readRDS("~/Dropbox/dds.RDS")
SummarizedExperiment::colData(dds)

normalized_counts <- DESeq2::normTransform(dds) %>% 
  SummarizedExperiment::assay()


colSums(normalized_counts) / 1e6

mod1 <- stats::model.matrix(~group_nr, data = metadata)

mod0 <- stats::model.matrix(~1, data = metadata)

svs <- sva::sva(normalized_counts, mod = mod1, mod0 = mod0)

svs$sv

cor(svs$sv[,1], mod1[,2])

DESeq2::plotPCA(DESeq2::normTransform(dds), intgroup = c("group_nr"))


RUVSeq::RUVg(normalized_counts, )


# PhD budget----
# Already planned + plane to Boston + Hotel in Whistler + Plane home + Hotel in Boston + Teaching Lab + Visualize your science + edx online course
25000+3000+14000+8000+7000+5500+5000+3000
