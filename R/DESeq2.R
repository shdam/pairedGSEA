#' Run DESeq2
#' 
#' @import DESeq2
#' @import sva
#' @export
runDESeq2 <- function(md,
                      archs4db,
                      groupCol,
                      comparison,
                      tpm = FALSE,
                      samples = "id",
                      prefilter = 10,
                      dds_out = FALSE,
                      parallel = FALSE,
                      cores = 4){
  
  # Look for database file
  if(!file.exists(archs4db)) stop("Database file is missing!\nLooking for: ", archs4db)
  # TPM file
  if(tpm & typeof(tpm) != "character"){
    archs4db_tpm <- archs4db %>% 
      stringr::str_replace("counts", "tpm")
  }
  if(tpm & !file.exists(archs4db_tpm)) stop("Database file is missing!\\nLooking for: ", archs4db_tpm)
  
  # Register parallel
  if(parallel & is.numeric(cores)) {
    missing_package("BiocParallel", "Bioc")
    BiocParallel::register(BiocParallel::MulticoreParam(4))
  }
  
  # Loading metadata
  metadata <- prepMeta(md_file, groupCol, comparison)
  
  # Define samples
  if(samples %in% colnames(metadata)) samples <- metadata[[samples]]
  else if(!(typeof(samples) == "character" & length(samples) > 1)) stop("Please specificy 'samples' as a column in metadata or as a vector of samples in database.")
  
  ### Load count matrix
  txCount <- loadArchs4(samples, archs4db)
  ### Ensure rows in metadata matches columns in the count matrix
  txCount <- txCount[, metadata$id]
  ### Pre-filter
  if(prefilter) txCount <- preFilter(txCount, prefilter)
  
  ### SVA
  dds <- runSVA(txCount, metadata, groupCol)
  
  message("Running DESeq2")
  
  dds <- DESeq2::DESeq(dds)
  
  if(dds_out) {
    check_make_dir("results/")
    saveRDS(dds, paste0("results/", dds_out))
  }
  # Ensure correct format for comparison
  if(typeof(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  # Extract results
  res <- DESeq2::results(dds, contrast = c(groupCol, comparison))
  
  # Add TPM
  if(tpm) res <- addTPM(res, samples, archs4db_tpm)
  
  return(res)
}
