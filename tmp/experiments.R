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
