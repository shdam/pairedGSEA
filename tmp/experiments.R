### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()



### Combine each and extract the comparisons to be run
experiments <- combineExperiments()

### Load GTF file
gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

BiocParallel::register(MulticoreParam(workers = 10))

# row <- experiments[82,]

### Run experiments ----

### Define file to read from and group column
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
groupCol <- "group_nr"

apply(experiments[1,], 1, runExperiment, archs4db)


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

