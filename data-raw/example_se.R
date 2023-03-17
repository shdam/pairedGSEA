## code to prepare `example_se` dataset goes here
pkgload::load_all(".")
library("magrittr")
library("dplyr")
library("stringr")

experiment_title <- "77_GSE61220_TNF Treatment 12hrs"

# Download results from Zenodo
url <- "https://zenodo.org/record/7032090/files/77_GSE61220_TNF%20Treatment%2012hrs.RDS?download=1"
download.file(url, destfile = "example_experiment.rds", method = "curl")

experiment <- readRDS("example_experiment.rds")
system2("rm", "example_experiment.rds")


# Identify "Telomer" gene sets and extract related genes
experiment_ora_small <- experiment$gene_sets %>% 
  filter(str_detect(tolower(pathway), tolower("Telomer")))

interesting_genes <- unlist(experiment_ora_small$overlapGenes_dexseq) %>% 
  unique()

# Add 900 random genes
set.seed(900)
random_genes <- experiment$genes %>% 
  slice_sample(n = 900) %>% 
  pull(gene) %>% 
  c(interesting_genes) %>% 
  unique()


# Prepare count matrix from ARCHS4 database
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

tx_count <- prepare_tx_count(
  metadata = experiment$metadata,
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
  colData = experiment$metadata %>% 
    dplyr::filter(id %in% colnames(tx_count_example)) %>% 
    dplyr::mutate(group_nr = factor(group_nr)) %>% 
    dplyr::select(study, id, source, final_description, group_nr) %>% 
    S4Vectors::DataFrame()
)

# Use data
usethis::use_data(example_se, overwrite = TRUE)
