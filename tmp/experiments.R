### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()



### Combine each and extract the comparisons to be run
experiments <- combineExperiments()

### Load GTF file
gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

BiocParallel::register(MulticoreParam(workers = 10))

# row <- experiments[19,]

### Run experiments ----

### Define file to read from and group column
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
groupCol <- "group_nr"

apply(experiments[82,], 1, runExperiment, archs4db)
row <- experiments[19,]

### Analyse experiment results ----

# Load MSigDB
gene_sets <- prepMsigdb()

apply(row, 1, analyseExperiment)
apply(experiments, 1, analyseExperiment)




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


concatResults <- function(experiments){
  
  concatFgsea <- tribble(~pathway, ~experiment, ~deseq2, ~dexseq, ~dexseq2, ~overlap)
  
  for(row in 1:nrow(experiments)){
    row <- experiments[row, ]
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
    fgseaRes <- readRDS(paste0("results/", dataname, "_fgseaRes_", experimentTitle, ".RDS"))
    fgseaDxr <- readRDS(paste0("results/", dataname, "_fgseaDxr_", experimentTitle, ".RDS"))
    fgseaDxr2 <- readRDS(paste0("results/", dataname, "_fgseaDxr2_", experimentTitle, ".RDS"))
    
    ### Significant gene sets
    pathRes <- fgseaRes %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathDxr <- fgseaDxr %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathDxr2 <- fgseaDxr2 %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    
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

concatForaResults <- function(experiments){
  
  concatFora <- tribble(~pathway, ~experiment, ~deseq2, ~dexseq, ~decombined)
  foratot <- tribble(~pathway, ~experiment, ~padj_deseq2, ~padj_dexseq, ~padj_decombined)
  
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
      rename(padj_deseq2 = padj) %>% 
      full_join(foradxr, by = "pathway") %>% 
      rename(padj_dexseq = padj) %>% 
      full_join(foraresdxr, by = "pathway") %>% 
      rename(padj_decombined = padj) %>% 
      select(pathway, starts_with("padj")) %>% 
      mutate(experiment = paste(row$study, experimentTitle))
    foratot <- foratot %>% 
      dplyr::bind_rows(comb)

    
    ### Significant gene sets
    pathres <- forares %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathdxr <- foradxr %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    pathresdxr <- foraresdxr %>% filter(padj < 0.05) %>% select(pathway) %>% as_tibble
    
    overlap <- union(pathres, pathdxr) %>% union(pathresdxr) %>% 
      mutate(deseq2 = pathway %in% pathres$pathway,
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

concatFora <- readRDS("results/concatFora.RDS")
foratot <- readRDS("results/foratot.RDS")

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

# Similarities
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

# Combine alle deseq2 and dexseq results
concatRes <- function(experiments){
  
  concatResults <- tibble::tribble(~gene, ~transcript, ~experiment, ~baseMean, ~log2FC, ~lfcSE, ~stat_deseq2,
                           ~pvalue_deseq2, ~padj_deseq2, ~tpm, ~exonBaseMean, ~dispersion, ~stat_dexseq, ~pvalue_dexseq,
                           ~padj_dexseq, ~log2FC_baseline_vs_condition)
  
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
      dplyr::select(colnames(concatResults))
    
    
    concatResults <- concatResults %>% 
      dplyr::bind_rows(comb)
  }
  saveRDS(concatResults, "results/concatResults.RDS")
  return(concatResults)
}
concatResults <- concatRes(experiments)

# combine all deseq2 results
concatRes <- function(experiments){
  
  concatResults <- tibble::tribble(~gene, ~transcript, ~experiment, ~baseMean, ~log2FC, ~lfcSE, ~stat,
                                   ~pvalue, ~padj, ~tpm)
  
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
      dplyr::select(colnames(concatResults))
    
    
    concatResults <- concatResults %>% 
      dplyr::bind_rows(comb)
  }
  saveRDS(concatResults, "results/concatResults.RDS")
  return(concatResults)
}

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
concatGene <- function(experiments){
  
  concatgene <- tibble::tibble()
  
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
    comb <- readRDS(paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
    comb$experiment <- paste(row$study, experimentTitle)
    
    
    concatgene <- concatgene %>% 
      dplyr::bind_rows(comb)
  }
  saveRDS(concatgene, "results/concatGenes.RDS")
  return(concatgene)
}

concatGenes <- concatGene(experiments)

umap_data <- concatGenes %>% 
  dplyr::select(experiment, ensembl_gene, lfc_deseq2) %>% 
  dplyr::rename(gene = ensembl_gene,
                lfc = lfc_deseq2) %>% 
  dplyr::filter(!is.na(lfc)) %>% 
  tidyr::pivot_wider(names_from = "experiment", values_from = "lfc")