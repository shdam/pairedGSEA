### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()

### List metadata files
md_files <- list.files("metadata", full.names = TRUE)

### Combine each and extract the comparisons to be run
experiments <- lapply(md_files, FUN = function(x) {df <- readxl::read_xlsx(x); df$filename <- x; df})  %>% 
  dplyr::bind_rows() %>% 
  dplyr::filter(!is.na(`comparison_title (empty_if_not_okay)`))




BiocParallel::register(MulticoreParam(workers = 10))

# row <- experiments[82,]

### Run experiments
runExperiment <- function(row){
  
  row <- tibble::as_tibble(row, rownames = "names") %>% 
    tidyr::pivot_wider(values_from = value, names_from = names)
  
  message("Running on ", row$study)
  
  ### Load metadata
  md_file <- row$filename
  dataname <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove("csv") 
  ### Define file to read from
  archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
  ### Define tpm file
  tpm <- stringr::str_replace(archs4db, "counts", "tpm")
  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)`
  experimentTitle <- row$`comparison_title (empty_if_not_okay)`
  groupCol <- "group_nr"
  
  
  ### Run in parallel
  # BiocParallel::register(BiocParallel::MulticoreParam(4))
  
  ### Prepare for DE
  dds <- prepDE(md = md_file,
                archs4db = archs4db,
                groupCol = groupCol,
                comparison = comparison,
                prefilter = 10)
  
  # return(dds)}
  
  ### Run DESeq2
  res_deseq2 <- runDESeq2(dds,
                          groupCol = groupCol,
                          comparison = comparison,
                          samples = dds$id,
                          tpm = tpm,
                          parallel = TRUE,
                          fitType = "local",
                          BPPARAM = BiocParallel::bpparam())#, dds_out = "deseq2_1_GSE154968.RDS")
  
  saveRDS(res_deseq2, paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
  rm(res_deseq2)
  
  ### Prepare DEXSeq
  dxr <- runDEXSeq(dds, groupCol, comparison)
  
  
  
  saveRDS(dxr, paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
  
  message(row$study, " is done.")
  
  
  # or to shrink log fold changes association with condition:
  #res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
  
  
}

apply(experiments[82,], 1, runExperiment)

