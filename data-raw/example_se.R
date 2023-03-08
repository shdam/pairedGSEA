## code to prepare `example_se` dataset goes here
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
source("tmp/run_experiment.R")
library("magrittr")
library("dplyr")
library("stringr")

example_experiment <- "77_GSE61220_TNF Treatment 12hrs"

# Extract experiment genes
concatenated_genes <- readRDS("results/concatenated_genes.RDS")
experiment_data <- concatenated_genes %>% 
  filter(experiment == example_experiment)
rm(concatenated_genes)

# Extract experiment gene sets
ora_all <- readRDS("results/ora_all.RDS")

experiment_ora <- ora_all %>% 
  filter(experiment == example_experiment)
rm(ora_all)

# Identify "Telomer" gene sets and extract related genes
experiment_ora_small <- experiment_ora %>% 
  filter(str_detect(tolower(pathway), tolower("Telomer")))

interesting_genes <- unlist(experiment_ora_small$overlapGenes_dexseq) %>% 
  unique()

# Add 900 random genes
set.seed(900)
random_genes <- experiment_data %>% 
  slice_sample(n = 900) %>% 
  pull(gene) %>% 
  c(interesting_genes) %>% 
  unique()


# Prepare count matrix
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

md_file <- "metadata/77_GSE61220.xlsx"
metadata <- readxl::read_excel(md_file)

tx_count <- prepare_tx_count(
  metadata = metadata,
  gtf = gtf,
  archs4db = archs4db,
  group_col = "group_nr",
  baseline_case = c(1, 2)
)

# Isolate rows related to genes of interest and random genes, remove those that will not pass pre-filtering step
tx_count_example <- tx_count[which((str_split(rownames(tx_count), pattern = ":", simplify = T)[, 1] %in% random_genes) &
                                     rowSums(tx_count) > 10),]

# Combine count matrix and metadata in Summarized Experiment
example_data <- SummarizedExperiment::SummarizedExperiment(
  assays = list("counts" = tx_count_example),
  colData = metadata %>% 
    dplyr::filter(id %in% colnames(tx_count_example)) %>% 
    dplyr::mutate(group_nr = factor(group_nr)) %>% 
    dplyr::select(study, id, source, final_description, group_nr) %>% 
    S4Vectors::DataFrame()
)

# Use data
usethis::use_data(example_se, overwrite = TRUE)
