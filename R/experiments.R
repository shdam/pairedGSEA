#!/home/ctools/opt/R-4.3.1/bin/Rscript
# run with: ./R/experiments.R -e 1:50 -w 5 -c 4
# Optparse and json libraries ----
library("optparse")
library("jsonlite")

# Script arguments
option_list <- list(
  optparse::make_option(
    c("-e", "--experiments"),
    type = "character", default = "all",
    help = "Experiments to run."
  ),
  optparse::make_option(
    c("--gsea"),
    action = "store_true",
    type = "logical", default = FALSE,
    help = "Only run GSEA."
  ),
  optparse::make_option(
    c("-g", "--gtf"),
    action = "store_true",
    type = "logical", default = FALSE,
    help = "Use GTF."
    ),
  optparse::make_option(
    c("-l", "--limma"),
    action = "store_true",
    type = "logical", default = FALSE,
    help = "Use limma."
  ),
  optparse::make_option(
    c("-s", "--sva"),
    action = "store_false",
    type = "logical", default = TRUE,
    help = "Run SVA"
  ),
  optparse::make_option(
    c("-v", "--verbose"),
    action = "store_true",
    type = "logical", default = FALSE,
    help = "Verbose."
  ),
  optparse::make_option(
    c("-w", "--workers"),
    type = "numeric", default = 6,
    help = "Workers to use in BPPARAM."
  ),
  optparse::make_option(
    c("-c", "--cores"),
    type = "numeric", default = 4,
    help = "Cores to use in pbmcmapply."
  ),
  optparse::make_option(
    c("-a", "--archs4version"),
    type = "character", default = "1.1",
    help = "ARCHS4 version to use."
  )
  )
opt_parser <- optparse::OptionParser(
  usage = "-e <e.g. 1:50> -g -w <workers> -c <cores>",
  option_list = option_list
)
opt <- optparse::parse_args(opt_parser)

# Define parameters ----
use_gtf <- opt$gtf
use_limma <- opt$limma
run_sva <- opt$sva
only_gsea <- opt$gsea
workers <- opt$workers
cores <- opt$cores
rows <- opt$experiments
verbose <- opt$verbose
archs4version <- opt$archs4version

if (interactive()) { # interactive
  rows <- 2
  use_gtf <- FALSE
  workers <- 4
  cores <- 4
  use_limma <- FALSE
  run_sva <- TRUE
  verbose <- TRUE
  quiet <- !verbose
  only_gsea <- FALSE
  archs4version <- "1.1"
}



### Load package ----
# pkgload::load_all(path = "/home/projects/shd/pairedGSEA")
# remotes::install_github("shdam/pairedGSEA", ref = "interactant_fcs")
# library("tidyverse")
lib.loc <- "/home/people/sohdam/R/x86_64-pc-linux-gnu-library/4.3/"
library("pairedGSEA", lib.loc = lib.loc)
library("magrittr")
library("ggplot2")
library("pbmcapply", lib.loc = lib.loc)



source("R/run_experiment.R")
source("R/run_analysis.R")
source("R/concat.R")

rm(aggregate_pvalue)
theme_set(theme_classic(base_size = 20))
### Combine each and extract the comparisons to be run

experiments <- combine_experiments(limma = use_limma)

if (rows == "all") 
  rows <- seq_len(nrow(experiments)) else rows <- eval(parse(text = rows))

### Load GTF file ----
if (use_gtf && !file.exists("gtfextract.rds")) {
  track <- rtracklayer::import("/home/databases/ensemble/Homo_sapiens.GRCh38.107.gtf.gz")
  gtf <- tibble::tibble(
    gene_id = track$gene_id, 
    transcript_id = track$transcript_id,
    gene_name = track$gene_name,
    gene_tx = stringr::str_c(gene_id, transcript_id, sep = ":")
  ) |> 
    unique()
  rm(track)
  saveRDS(gtf, "gtfextract.rds")
} else if (use_gtf){
  gtf <- readRDS("gtfextract.rds")
  # gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)
} else {
  gtf <- NULL
}

### Define file to read from
if (archs4version == "1.1") {
  archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
} else if (archs4version == "2.2") {
  archs4db <- "/home/databases/archs4/v2.2/human_transcript_v2.2.h5"
} else {
  stop("Only versions '2.2' and '1.1' are usable. Default is '1.1'.")
}

### Register cores
BiocParallel::register(BiocParallel::MulticoreParam(workers = workers))

### Run experiments ----
if (FALSE){
  #tests
  row <- experiments[1,]
  apply(experiments[1, ], 1, run_experiment, archs4db)
  
  gene_sets <- pairedGSEA::prepare_msigdb()
  apply(experiments[1, ], 1, run_analysis, gene_sets, run_fgsea = TRUE)
}

if (!only_gsea) {
  message("Running ", length(rows), " analyses using ", workers, " workers and ", cores, " cores", ":")
  
  empty <- pbmcapply::pbmclapply(
    rows, run_experiment, experiments = experiments, archs4db = archs4db, quiet = !verbose,
    mc.cores = cores, ignore.interactive = TRUE
  )
  print(empty)
  stopifnot("Not all analyses could be completed" = length(empty) == length(rows))
  
  message("Paired differential analyses are done!")
}

# Run GSEA ----
message("Initiating Gene-Set Enrichment Analyses..")

gene_sets <- pairedGSEA::prepare_msigdb()
empty <- pbmcapply::pbmclapply(
  rows, run_analysis, experiments = experiments, gene_sets = gene_sets, run_fgsea = TRUE, quiet = !verbose,
  mc.cores = cores, ignore.interactive = TRUE
  )
print(empty)
message("Done!")
