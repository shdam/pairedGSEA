### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()

source("tmp/run_experiment.R")
source("tmp/run_analysis.R")
source("tmp/concat.R")

rm(aggregate_pvalue)
### Combine each and extract the comparisons to be run
experiments <- combine_experiments()

### Load GTF file
# gtf <- readRDS("gtfextract.rds")
# gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

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

apply(experiments[1, ], 1, run_experiment, archs4db, gtf = NULL)
apply(experiments, 1, run_analysis, gene_sets)
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


# Concat gene level ----

# concatenated_genes <- concatenate_genes(experiments)

concatenated_genes <- readRDS("results/concatenated_genes.RDS")

# Median genes per analysis ----
concatenated_genes %>% 
  group_by(experiment) %>% 
  summarise(deseq = sum(padj_deseq < 0.05, na.rm = TRUE), dexseq = sum(padj_dexseq < 0.05, na.rm = TRUE)) %>% 
  summarise(med_deseq = median(deseq), med_dexseq = median(dexseq))




umap_data <- concatenated_genes %>% 
  dplyr::select(experiment, ensembl_gene, lfc_deseq2) %>% 
  dplyr::rename(gene = ensembl_gene,
                lfc = lfc_deseq2) %>% 
  dplyr::filter(!is.na(lfc)) %>% 
  dplyr::distinct(gene, experiment, .keep_all = TRUE) %>% 
  tidyr::pivot_wider(names_from = "experiment", values_from = "lfc")

um <- umap_data %>% 
  dplyr::select(-gene) %>% 
  dplyr::filter(across(everything(), ~ !is.na(.x))) %>% 
  as.matrix() %>% 
  t() %>% 
  uwot::umap(n_neighbors = 3)

p3 <- um %>% 
  as_tibble(rownames = "Experiment") %>% 
  ggplot(aes(x = V1, y = V2))+
  geom_point() + 
  labs(caption = "n_neighbors = 3")



# SVA metrics -------------------------------------------------------------

# concatenated_design <- concatenate_design(experiments)

(total_svas <- concatenated_design %>% 
  dplyr::mutate(design = as.character(design),
                svs = stringr::str_count(design, "sv\\d")) %>% 
  dplyr::count(svs))



# Fig 1b ------------------------------------------------------------------

library("dplyr")
library("tidyr")
library("ggplot2")
concatenated_genes <- readRDS("results/concatenated_genes.RDS")
concatenated_genes <- readRDS("results/concatGenes_small.RDS")

# concatenated_genes <- concatenated_genes %>% 
#     filter(pvalue_deseq < 0.05 | pvalue_dexseq < 0.05)
# saveRDS(concatenated_genes, "results/cconcatenated_genes_small.RDS")

# Number of genes

found_genes <- concatenated_genes %>% 
  group_by(experiment) %>%
  summarise(`Differential\nExpression` = sum(padj_deseq < 0.05, na.rm = TRUE),
            `Differential\nSplicing` = sum(padj_dexseq < 0.05, na.rm = TRUE),
            Overlap = sum(padj_dexseq < 0.05 & padj_deseq < 0.05, na.rm = TRUE),
            `Only\nDifferential\nSplicing` = sum((padj_dexseq < 0.05) & !(padj_deseq < 0.05), na.rm = TRUE))

## as density
found_genes %>% 
  pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "genes") %>% 
  ggplot(aes(x = genes, fill = Analysis)) +
  geom_density(alpha = 0.5)

## as violin
fig1b <- found_genes %>% 
  pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Genes") %>% 
  mutate(Analysis = factor(Analysis, levels = c("Differential\nExpression",
                                                "Differential\nSplicing",
                                                "Overlap",
                                                "Only\nDifferential\nSplicing"))) %>% 
  ggplot(aes(y = Genes, x = Analysis)) +
  scale_y_log10() +
  ggforce::geom_sina() +
  # stat_summary(fun = median, geom = "point", size = 3, color = "red", shape = 18) +
  # geom_violin() +
  geom_boxplot(width = 0.05) +
  labs(x = "",
       y = "Significant genes per experiment") +
  theme_classic(base_size = 20)

ggsave(plot = fig1b, filename = "figs/fig1b.png")


# Fig 1c ------------------------------------------------------------------

gene_similarity <- concatenated_genes %>%
  mutate(DESeq2 = padj_deseq < 0.05,
         DEXSeq = padj_dexseq < 0.05) %>% 
  group_by(experiment) %>% 
  group_map(.f = function(df, ...){
    simpson <- df %>% 
      select(starts_with("DE")) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    # jaccard <- df %>% 
    #   select(starts_with("DE")) %>% 
    #   proxy::simil(by_rows = FALSE, method = "Jaccard")
    tibble(`Overlap Coefficient` = simpson[1]#, Jaccard = jaccard[1]
           )
  }) %>% 
  bind_rows() %>% 
  mutate(experiment = concatenated_genes %>% distinct(experiment) %>% pull())

fig1c <- gene_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  ggplot(aes(x = Similarity, fill = Method)) +
  geom_density(alpha = 0.5) +
  # geom_histogram(binwidth = 0.02) +
  labs(#title = "Overlap coefficient between DGE and DGS significant genes",
       x = "Overlap coefficient",
       y = "Density") +
  scale_fill_manual(values = "gray") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

ggsave(plot = fig1c, filename = "figs/fig1c.png")

# Fig 1d ------------------------------------------------------------------

# Distribution of fraction of tested transcripts that are DTE? 
# (with the notion if it is only 1 of 10 transcripts that are DTE 
# it is not the gene that is differentially expressed but a single isoofrm.
# Vise verca if all isoforms are upregulated it is DGE not DGU).

concatenated_genes <- readRDS("results/concatenated_genes.RDS")

significant_genes <- concatenated_genes %>% 
  filter(padj_dexseq < 0.05 & padj_deseq < 0.05)

# concatenated_results <- concatenate_results(experiments)
concatenated_results <- readRDS("results/concatenated_results.RDS")

# analysed_genes <- concatenated_results %>% 
#   count(experiment, gene) %>% 
#   filter(n > 1) %>% 
#   distinct(experiment, gene)

transcript_fractions <- concatenated_results %>% 
  semi_join(significant_genes, by = c("experiment", "gene")) %>% 
  # filter(gene %in% analysed_genes) %>% 
  group_by(experiment, gene) %>% 
  summarise(fraction = sum(padj_dexseq < 0.05, na.rm = TRUE) / n(), n = n()) %>% 
  filter(n > 1)

saveRDS(transcript_fractions, "results/transcript_fractions.RDS")

transcript_fractions <- readRDS("results/transcript_fractions.RDS")

exclude <- transcript_fractions %>% 
  count(experiment) %>% 
  filter(n < 10)

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

fig1d <- transcript_fractions %>% 
  # anti_join(exclude, by = "experiment") %>% 
  filter(gene != "NA") %>% 
  ggplot() +
  aes(x = fraction*n, y = n, color = experiment) +
  geom_point(alpha = 0.1) +
  # geom_density(fill = NA) +
  scale_color_manual(values = rep("gray", length(unique(transcript_fractions$experiment)))) +
  labs(y = "Number of transcripts in gene",
       x = "Fraction of significant transcripts in a gene") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

# Fig 2a ------------------------------------------------------------------


ora_all <- readRDS("results/ora_all.RDS")
ora_all <- rename(ora_all, experiment = experiment_title)

# ora_all <- ora_all %>% 
#   select(-starts_with("overlapG")) %>% 
#   filter(padj_deseq < 0.05 | padj_dexseq < 0.05 | padj_decombined < 0.05)
# saveRDS(ora_all, "results/ora_all_small.RDS")
# Number of significant pathways

found_pathways <- ora_all %>% 
  group_by(experiment) %>% 
  summarise(`Differential\nExpression` = sum(padj_deseq < 0.05, na.rm = TRUE),
            `Differential\nSplicing` = sum(padj_dexseq < 0.05, na.rm = TRUE),
            Overlap = sum((padj_dexseq < 0.05) & (padj_deseq < 0.05), na.rm = TRUE),
            # Paired = sum(padj_decombined < 0.05, na.rm = TRUE),
            # Overlap2 = sum((padj_decombined < 0.05) & (padj_deseq2 < 0.05), na.rm = TRUE),
            # `Only paired` = sum((padj_decombined < 0.05) & !(padj_deseq2 < 0.05), na.rm = TRUE),
            # `Only DESeq2` = sum(!(padj_decombined < 0.05) & (padj_deseq2 < 0.05), na.rm = TRUE),
            `Only\nDifferential\nSplicing` = sum(!(padj_deseq < 0.05) & (padj_dexseq < 0.05), na.rm = TRUE)
            )

## as violin
fig2a <- found_pathways %>% 
  pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Pathways") %>% 
  # filter(Pathways != 0) %>%
  mutate(Analysis = factor(Analysis, levels = c("Differential\nExpression",
                                                "Differential\nSplicing",
                                                "Overlap",
                                                "Only\nDifferential\nSplicing")
                             #c("DESeq2", "Paired", "Overlap", "Only paired", "Only DESeq2")
                             )) %>% 
  ggplot(aes(y = Pathways, x = Analysis)) +
  scale_y_log10() +
  # stat_summary(fun = median, geom = "point", size = 3, color = "red", shape = 18) +
  # geom_violin() +
  ggforce::geom_sina() +
  geom_boxplot(width = 0.05) +
  labs(x = "",
       y = "Significant pathways per experiment") +
  theme_classic(base_size = 20)

ggsave(plot = fig2a, filename = "figs/fig2a.png")

# Fig 2b ------------------------------------------------------------------


pathway_similarity <- ora_all %>%
  mutate(DESeq2 = padj_deseq < 0.05,
         DEXSeq = padj_dexseq < 0.05,
         # paired = padj_paired < 0.05,
         added = padj_deseq < 0.05 | padj_dexseq < 0.05) %>% 
  group_by(experiment) %>% 
  group_map(.f = function(df, ...){
    dge_dtu <- df %>% 
      select(DESeq2, DEXSeq) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    # dge_paired <- df %>% 
    #   select(DESeq2, paired) %>% 
    #   proxy::simil(by_rows = FALSE, method = "Simpson")
    # paired_added <- df %>% 
    #   select(added, paired) %>% 
    #   proxy::simil(by_rows = FALSE, method = "Simpson")
    tibble(DGE_DTU = dge_dtu[1]
           # , DGE_Paired = dge_paired[1], Paired_Added = paired_added[1]
           )
  }) %>% 
  bind_rows() %>% 
  mutate(experiment = ora_all %>% distinct(experiment) %>% pull())


# Median line
eta <- pathway_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  filter(!is.na(Similarity)) %>% 
  group_by(Method) %>% 
  summarise(group_median = median(Similarity))

# as density
fig2b <- pathway_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  ggplot(aes(x = Similarity, fill = Method)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = eta, aes(xintercept = group_median, color = Method),
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

# Fig 2e ------------------------------------------------------------------

# ora_all <- concatenate_ora(experiments)
ora_all <- readRDS("results/ora_all.RDS")

ora_all %>% 
  filter(padj_dexseq < 0.05) %>%
  ggplot(aes(x = relative_risk_dexseq,
             y = relative_risk_deseq)) +
  geom_point()

ora_all %>% 
  filter(padj_dexseq < 0.05,
         abs(relative_risk_dexseq - relative_risk_deseq) > 1) %>%
  ggplot(aes(x = relative_risk_dexseq - relative_risk_deseq)) +
  geom_density()

fig2e <- ora_all %>% 
  filter(padj_dexseq < 0.05) %>%
  mutate(enrichment_score_s = relative_risk_dexseq - relative_risk_deseq) %>% 
  filter(abs(enrichment_score_s) > 1) %>% 
  mutate(enrichment_score_s = log2(abs(enrichment_score_s))* sign(enrichment_score_s)) %>%
  ggplot() +
  aes(x = enrichment_score_s, fill = experiment) +
  # geom_density(alpha = 0.02, color = NA) +
  geom_density(fill = NA) +
  scale_fill_manual(values = rep("gray", length(unique(ora_all$experiment)))) +
  # scale_x_continuous(trans = "log2") +
  # coord_cartesian(xlim = c(-4,10)) +
  labs(y = "Density",
       x = "Enrichment Score Shift") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

ggsave(plot = fig2e, filename = "figs/fig2e.png")

# Pathway overlap ----

ora_all %>% 
  mutate(DESeq2 = overlap_deseq / size,
         DEXSeq = overlap_dexseq / size#,
         # Combined = overlap_decombined / size
         ) %>% 
  select(pathway, DESeq2, DEXSeq
         # , Combined
         ) %>% 
  pivot_longer(cols = -pathway, names_to = "Method", values_to = "Overlap") %>% 
  filter(!is.na(Overlap)) %>% 
  ggplot(aes(x = Overlap, color = Method)) +
  geom_density(alpha = 0) +
  labs(x = "Pathway overlap (%)", y = "Density") +
  theme_classic()


# Odds ratio ----


ora_res <- readRDS("results/1_GSE154968_HaCaT cells transfected with miR-21-3p_ora_deseq.RDS")
overlap <- 277
size <- 990
size_genes <- 1748
size_universe <- 17073

(overlap / size_genes) / (size / size_universe)
# (overlap / size) / (size_genes / size_universe)
# (overlap*size_universe)/(size_genes*size)
ora_res %>% mutate(odds_ratio = (overlap / size_genes) / (size_geneset / size_universe)) %>% pull(odds_ratio) %>% qplot()
# ora_res %>% mutate(odds_ratio = (overlap / size_geneset) / (size_genes / size_universe)) %>% pull(odds_ratio) %>% qplot()
## Gene correlation ----

correlated_genes <- concatGenes %>% 
  filter(pvalue_deseq2 < 0.05 | pvalue_dexseq < 0.05,
         !is.na(pvalue_deseq2),
         !is.na(pvalue_dexseq)) %>% 
  group_by(experiment) %>% 
  # mutate(sig_deseq2 = pvalue_deseq2 < 0.05,
  #        sig_dexseq = pvalue_dexseq < 0.05) %>% 
  summarise(correlation = cor(pvalue_deseq2, pvalue_dexseq))

qplot(correlated_genes$correlation)

## fgsa analysis ----

concatFgsea <- concatFgseaResults(experiments)

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

# V1
rr_shifts %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(DESeq2 = relative_risk_deseq,
                DEXSeq = relative_risk_dexseq,
                Shift = relative_risk_shift) %>% 
  pivot_longer(cols = -c("experiment", "pathway"), names_to = "Analysis", values_to = "Pathways") %>% 
  mutate(Analysis = factor(Analysis, levels = c("DESeq2", "DEXSeq", "Shift")
                           #c("DESeq2", "Paired", "Overlap", "Only paired", "Only DESeq2")
  )) %>% 
  ggplot(aes(y = Pathways, x = Analysis)) +
  # stat_summary(fun = median, geom = "point", size = 3, color = "red", shape = 18) +
  # geom_violin() +
  ggforce::geom_sina() +
  # geom_boxplot(width = 0.05) +
  labs(x = "",
       y = "Odds Ratio") +
  theme_classic()

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
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low = "darkred", high = "navy")


experiment_title <- "76_GSE183984_Cetuximab treatment of primary colorectal cancer"
# experiment_title <- "79_GSE139262_SMARCB1 overexpression"

dexseqres <- readRDS(paste0("results/", experiment_title, "_dexseqres.RDS"))
deseqres <- readRDS(paste0("results/", experiment_title, "_deseq2res.RDS"))
ora_dexseq <- readRDS(paste0("results/", experiment_title, "_ora_dexseq.RDS"))
ora_deseq <- readRDS(paste0("results/", experiment_title, "_ora_deseq.RDS"))
aggregated_pvals <- readRDS(paste0("results/", experiment_title, "_aggregated_pvals.RDS"))

# PhD budget----
# Already planned + plane to Boston + Hotel in Whistler + Plane home + Hotel in Boston + Teaching Lab + Visualize your science + edx online course
25000+3000+14000+8000+7000+5500+5000+3000
