### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()



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
gene_sets <- prepare_msigdb()

apply(row, 1, run_experiment, archs4db)
apply(experiments, 1, analyseExperiment)
apply(row, 1, analyse_experiment)
apply(experiments, 1, getDDS, archs4db)



### Check missing
test <- stringr::str_c(experiments$study, experiments$`comparison_title (empty_if_not_okay)`)

test[which(table(test)>1)]

apply(experiments[rerun, ], 1, runExperiment)
apply(experiments[rerun, ], 1, analyseExperiment)

named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

rerun <- c(59,60, 73, 74, 184, 185)

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



### Fora analysis ----



concatFora <- concatForaResults(experiments)
concatFora <- readRDS("results/concatFora.RDS")
foratot <- readRDS("results/foratot.RDS")
foratot <- foratot %>% 
  dplyr::select(-dplyr::starts_with("overlapG"))
concatFora %>% 
  # filter(experiment == experiment[[1]]) %>% 
  pivot_longer(cols = starts_with("de"),names_to = "analysis", values_to = "present") %>% 
  mutate(present = as.factor(present),
         id = 1:nrow(.)) %>% 
  ggplot(aes(x = analysis, y = id, fill = present)) +
  geom_tile()


concatFora %>% 
  filter(decombined & !dexseq & !deseq2) %>% 
  count(experiment)


Study <- "GSE101968 Obstructed defecation symdrome"
concatFora %>% 
  filter(experiment == Study,
         decombined & !dexseq & !deseq2)

concatFora %>% 
  filter(!decombined & dexseq & !deseq2) %>% 
  count(experiment)
concatFora %>% 
  filter(!decombined & !dexseq & deseq2) %>% 
  count(experiment)

# Fora Similarities ----
concatFora %>% 
  group_by(experiment) %>% 
  summarise(jacdedex = proxy::simil(deseq2, dexseq, by_rows = FALSE, pairwise = TRUE, method = "Jaccard"),
            jacdecom = proxy::simil(deseq2, decombined, by_rows = FALSE, pairwise = TRUE, method = "Jaccard"),
            jacdexcom = proxy::simil(dexseq, decombined, by_rows = FALSE, pairwise = TRUE, method = "Jaccard"),
            simdedex = proxy::simil(deseq2, dexseq, by_rows = FALSE, pairwise = TRUE, method = "Simpson"),
            simdecom = proxy::simil(deseq2, decombined, by_rows = FALSE, pairwise = TRUE, method = "Simpson"),
            simdexcom = proxy::simil(dexseq, decombined, by_rows = FALSE, pairwise = TRUE, method = "Simpson")
            )

# Pathways per analysis
concatFora %>% 
  group_by(experiment) %>% 
  summarise(deseq2 = sum(deseq2), dexseq = sum(dexseq), decombined = sum(decombined))

# Median pathways per analysis
concatFora %>% 
  group_by(experiment) %>% 
  summarise(deseq2 = sum(deseq2), dexseq = sum(dexseq), decombined = sum(decombined)) %>% 
  summarise(med_deseq2 = median(deseq2), med_dexseq = median(dexseq), med_decombined = median(decombined))

# Median combined not in other
concatFora %>% 
  filter((decombined & !dexseq & !deseq2) | (!decombined & (dexseq | dexseq))) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(decombined)) %>% 
  summarise(median = median(n))

# Median DEXSeq not in other
concatFora %>% 
  filter((dexseq & !deseq2) | (!dexseq & !dexseq)) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(dexseq)) %>% 
  summarise(median = median(n))

# Median deseq2 not in other
concatFora %>% 
  filter((!dexseq & deseq2) | (!dexseq & !dexseq)) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(deseq2)) %>% 
  summarise(median = median(n))

# Antallet af gene-sets hvor p-værdien er mindre end for hver af de individuelle analyser
foratot %>% 
  filter(across(starts_with("padj"), ~ .x < 0.05)) %>% 
  group_by(experiment) %>% 
  summarise(n = sum(padj_decombined < padj_dexseq & padj_decombined < padj_deseq2))


## Rank shifts ----

# Plot distribution of adjusted p-values
foratot %>% 
  dplyr::filter(padj_deseq2 < 0.05 & padj_decombined < 0.05) %>% 
  tidyr::pivot_longer(cols = c("padj_deseq2", "padj_decombined"), names_to = "analysis", values_to = "padj") %>% 
  ggplot(aes(x = padj, fill = analysis))+
  geom_density(alpha = 0.5)

# Plot distribution of coverage (overlap / size)
foratot %>% 
  dplyr::filter(padj_deseq2 < 0.05 | padj_decombined < 0.05) %>% 
  tidyr::pivot_longer(cols = c("overlap_deseq2", "overlap_decombined"), names_to = "analysis", values_to = "overlap") %>% 
  dplyr::mutate(coverage = overlap / size) %>% 
  ggplot(aes(x = coverage, fill = analysis))+
  geom_density(alpha = 0.5)


# > foratot %>% filter(padj_decombined < 0.05) %>% nrow()
# [1] 377792
# > foratot %>% filter(padj_deseq2 < 0.05) %>% nrow()
# [1] 370800
# > foratot %>% filter(padj_decombined < 0.1) %>% nrow()
# [1] 462871
# > foratot %>% filter(padj_deseq2 < 0.1) %>% nrow()
# [1] 451905
# > foratot %>% filter(padj_deseq2 == 1) %>% nrow()
# [1] 120395
# > foratot %>% filter(padj_decombined == 1) %>% nrow()
# [1] 71002

# Rank each geneset within the experiments
ranks <- foratot %>% 
  dplyr::filter(padj_deseq2 < 0.05 | padj_decombined < 0.05) %>% 
  dplyr::group_by(experiment) %>% 
  dplyr::mutate(rank_deseq2 = rank(pval_deseq2),
                rank_dexseq = rank(pval_dexseq),
                rank_combined = rank(pval_decombined),
                rank_shift = rank_deseq2 - rank_combined,
                rank_shift2 = rank_deseq2 - rank_dexseq) %>% 
  dplyr::arrange(desc(rank_shift))
qplot(ranks$rank_shift)

# Store a list of all the top 20 ranked pathways of the paired analysis that were not in the top 20 for the deseq2 analysis
ranks %>% 
  dplyr::filter(
    rank_deseq2 > 50 & rank_combined < 50
  ) %>%  dplyr::select(pathway, experiment, dplyr::starts_with("rank")) %>% dplyr::arrange(experiment) %>% 
  # pull(rank_shift) %>% qplot()
  write.csv("results/fora_combined_top50.csv")

# List of all the top 20 ranked pathways of the deseq2 analysis that were not in the top 20 for the paired analysis
ranks %>% 
  dplyr::filter(
    rank_deseq2 < 20 & rank_combined > 20
  ) %>%  dplyr::select(pathway, experiment, dplyr::starts_with("rank")) %>% dplyr::arrange(experiment)

ranks %>% 
  dplyr::filter(
    rank_deseq2 > 50 & rank_dexseq < 50
  ) %>%  dplyr::select(pathway, experiment, dplyr::starts_with("rank")) %>% 
  dplyr::arrange(dplyr::desc(rank_shift2)) %>% 
  # pull(rank_shift) %>% qplot()
  write.csv("results/fora_compared_top50.csv")

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

# concatGenes <- concatGene(experiments)

concatGenes <- readRDS("results/concatGenes.RDS")

umap_data <- concatGenes %>% 
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

# concatsva <- concatSVA(experiments)

(total_svas <- concatsva %>% 
  dplyr::mutate(design = as.character(design),
                svs = stringr::str_count(design, "sv")) %>% 
  dplyr::count(svs))



# Fig 1b ------------------------------------------------------------------

library("dplyr")
library("tidyr")
library("ggplot2")
concatGenes <- readRDS("results/concatGenes.RDS")
concatGenes <- readRDS("results/concatGenes_small.RDS")

# concatGenes <- concatGenes %>% 
#     filter(pvalue_deseq2 < 0.05 | pvalue_dexseq < 0.05)
# saveRDS(concatGenes, "results/concatGenes_small.RDS")

# Number of genes
concatGenes <- concatGenes %>% 
  group_by(experiment) %>% 
  mutate(padj_deseq2 = p.adjust(pvalue_deseq2, "fdr"),
         padj_dexseq = p.adjust(pvalue_dexseq, "fdr"))
found_genes <- concatGenes %>% 
  # group_by(experiment) %>% 
  summarise(DESeq2 = sum(padj_deseq2 < 0.05, na.rm = TRUE),
            DEXSeq = sum(padj_dexseq < 0.05, na.rm = TRUE),
            Overlap = sum(padj_dexseq < 0.05 & padj_deseq2 < 0.05, na.rm = TRUE),
            `DEXSeq Only` = sum((padj_dexseq < 0.05) & !(padj_deseq2 < 0.05), na.rm = TRUE))

## as density
found_genes %>% 
  pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "genes") %>% 
  ggplot(aes(x = genes, fill = Analysis)) +
  geom_density(alpha = 0.5)

## as violin
fig1b <- found_genes %>% 
  pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Genes") %>% 
  mutate(Analysis = factor(Analysis, levels = c("DESeq2", "DEXSeq", "Overlap", "DEXSeq Only"))) %>% 
  ggplot(aes(y = Genes, x = Analysis)) +
  scale_y_log10() +
  ggforce::geom_sina() +
  # stat_summary(fun = median, geom = "point", size = 3, color = "red", shape = 18) +
  # geom_violin() +
  geom_boxplot(width = 0.05) +
  labs(x = "",
       y = "Significant genes per experiment") +
  theme_classic()

ggsave(plot = fig1b, filename = "figs/fig1b.png")


# Fig 1c ------------------------------------------------------------------

gene_similarity <- concatGenes %>%
  mutate(DESeq2 = padj_deseq2 < 0.05,
         DEXSeq = padj_dexseq < 0.05) %>% 
  group_by(experiment) %>% 
  group_map(.f = function(df, ...){
    simpson <- df %>% 
      select(starts_with("DE")) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    jaccard <- df %>% 
      select(starts_with("DE")) %>% 
      proxy::simil(by_rows = FALSE, method = "Jaccard")
    tibble(Simpson = simpson[1], Jaccard = jaccard[1])
  }) %>% 
  bind_rows() %>% 
  mutate(experiment = concatGenes %>% distinct(experiment) %>% pull())

fig1c <- gene_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  ggplot(aes(x = Similarity, fill = Method)) +
  geom_density(alpha = 0.5) +
  labs(x = "Similarity scores",
       y = "Density") +
  theme_classic()

ggsave(plot = fig1c, filename = "figs/fig1c.png")


# Fig 1d ------------------------------------------------------------------


foratot <- readRDS("results/foratot.RDS")

# foratot <- foratot %>% 
#   select(-starts_with("overlapG")) %>% 
#   filter(padj_deseq2 < 0.05 | padj_dexseq < 0.05 | padj_decombined < 0.05)
# saveRDS(foratot, "results/foratot_small.RDS")
# Number of significant pathways

found_pathways <- foratot %>% 
  group_by(experiment) %>% 
  summarise(DESeq2 = sum(padj_deseq2 < 0.05, na.rm = TRUE),
            # DEXSeq = sum(padj_dexseq < 0.05, na.rm = TRUE),
            Paired = sum(padj_decombined < 0.05, na.rm = TRUE),
            Overlap = sum((padj_decombined < 0.05) & (padj_deseq2 < 0.05), na.rm = TRUE),
            `Only paired` = sum((padj_decombined < 0.05) & !(padj_deseq2 < 0.05), na.rm = TRUE),
            `Only DESeq2` = sum(!(padj_decombined < 0.05) & (padj_deseq2 < 0.05), na.rm = TRUE)
            )

## as violin
fig1d <- found_pathways %>% 
  pivot_longer(cols = -experiment, names_to = "Analysis", values_to = "Pathways") %>% 
  # filter(Pathways != 0) %>%
  mutate(Analysis = factor(Analysis, levels = c("DESeq2", "Paired", "Overlap", "Only paired", "Only DESeq2"))) %>% 
  ggplot(aes(y = Pathways, x = Analysis)) +
  scale_y_log10() +
  # stat_summary(fun = median, geom = "point", size = 3, color = "red", shape = 18) +
  # geom_violin() +
  ggforce::geom_sina() +
  geom_boxplot(width = 0.05) +
  labs(x = "",
       y = "Significant pathways per experiment") +
  theme_classic()

ggsave(plot = fig1d, filename = "figs/fig1d.png")

# Fig 1e ------------------------------------------------------------------


pathway_similarity <- foratot %>%
  mutate(DESeq2 = padj_deseq2 < 0.05,
         DEXSeq = padj_dexseq < 0.05,
         paired = padj_decombined < 0.05,
         added = padj_deseq2 < 0.05 | padj_dexseq < 0.05) %>% 
  group_by(experiment) %>% 
  group_map(.f = function(df, ...){
    dge_dtu <- df %>% 
      select(DESeq2, DEXSeq) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    dge_paired <- df %>% 
      select(DESeq2, paired) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    paired_added <- df %>% 
      select(added, paired) %>% 
      proxy::simil(by_rows = FALSE, method = "Simpson")
    tibble(DGE_DTU = dge_dtu[1], DGE_Paired = dge_paired[1], Paired_Added = paired_added[1])
  }) %>% 
  bind_rows() %>% 
  mutate(experiment = foratot %>% distinct(experiment) %>% pull())


# Median line
eta <- pathway_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  filter(!is.na(Similarity)) %>% 
  group_by(Method) %>% 
  summarise(group_median = median(Similarity))

# as density
fig1e <- pathway_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  ggplot(aes(x = Similarity, color = Method)) +
  geom_density(alpha = 0) +
  geom_vline(data = eta, aes(xintercept = group_median, color = Method),
             linetype = "dashed") +
  # scale_color_brewer(palette = "BuPu") +
  labs(x = "Simpson Similarity") +
  theme_classic()

ggsave(plot = fig1e, filename = "figs/fig1e.png")

# as sina
pathway_similarity %>% 
  pivot_longer(cols = -experiment, names_to = "Method", values_to = "Similarity") %>% 
  filter(!is.na("Similarity")) %>% 
  ggplot(aes(y = Similarity, x = Method)) +
  ggforce::geom_sina() +
  theme_classic()

# Pathway overlap ----

foratot %>% 
  mutate(DESeq2 = overlap_deseq2 / size,
         DEXSeq = overlap_dexseq / size,
         Combined = overlap_decombined / size) %>% 
  select(pathway, DESeq2, DEXSeq, Combined) %>% 
  pivot_longer(cols = -pathway, names_to = "Method", values_to = "Overlap") %>% 
  filter(!is.na(Overlap)) %>% 
  ggplot(aes(x = Overlap, color = Method)) +
  geom_density(alpha = 0) +
  labs(x = "Pathway overlap (%)", y = "Density") +
  theme_classic()


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

# PhD budget
# Already planned + plane to Boston + Hotel in Whistler + Plane home + Hotel in Boston + Teaching Lab + Visualize your science + edx online course
25000+3000+14000+8000+7000+5500+5000+3000
