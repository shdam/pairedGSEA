## code to prepare `example` dataset goes here
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
source("tmp/run_experiment.R")
library(magrittr)
library(dplyr)
library(stringr)

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


experiment_ora_small <- experiment_ora %>% 
  filter(str_detect(tolower(pathway), tolower("Telomer")) |
         str_detect(tolower(pathway), tolower("Response")) |
         str_detect(tolower(pathway), tolower("POSITIVE")) |
         str_detect(tolower(pathway), tolower("GOBP_REGULATION")) |
         str_detect(tolower(pathway), tolower("GOBP_DNA_REPAIR")) |
         str_detect(tolower(pathway), tolower("GOCC_NUCLEAR_PROT")) |
         str_detect(tolower(pathway), tolower("GOBP_RNA_PROCESSI")) |
         str_detect(tolower(pathway), tolower("Hydro")))

interesting_genes <- unlist(experiment_ora_small$overlapGenes_dexseq) %>% 
  unique()

experiment_data_small <- experiment_data %>% 
  filter(gene %in% interesting_genes)



# Prepare count matrix



archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)




md_file <- "metadata/77_GSE61220.xlsx"

tx_count <- prepare_tx_count(
  metadata = md_file,
  gtf = gtf,
  archs4db = archs4db,
  group_col = "group_nr",
  baseline_case = c(1, 2)
)

tx_count_example <- tx_count[which(str_split(rownames(tx_count), pattern = ":", simplify = T)[, 1] %in% interesting_genes), ]


example_data <- list(
  "tx_count" = tx_count_example,
  "metadata" = readxl::read_excel(md_file)
)

usethis::use_data(example_data, overwrite = TRUE)
