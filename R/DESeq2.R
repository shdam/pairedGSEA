#' Pepare DE analysis
#' 
#' @import DESeq2
#' @import sva
#' @export
prepDE <- function(md,
                   archs4db,
                   groupCol,
                   comparison,
                   tpm = FALSE,
                   samples = "id",
                   prefilter = 10,
                   parallel = FALSE,
                   cores = 4){
  
  # Look for database file
  if(!file.exists(archs4db)) stop("Database file is missing!\nLooking for: ", archs4db)
  
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
 
  return(dds) 
}
#' Run DESeq2
#' 
#' @import DESeq2
#' @import sva
#' @export
runDESeq2 <- function(dds,
                      groupCol,
                      comparison,
                      tpm = FALSE,
                      samples = NA,
                      dds_out = FALSE,
                      parallel = FALSE,
                      cores = 4){
  # Register parallel
  if(parallel & is.numeric(cores)) {
    missing_package("BiocParallel", "Bioc")
    BiocParallel::register(BiocParallel::MulticoreParam(4))
  }
  message("Running DESeq2")
  dds <- DESeq2::DESeq(dds)
  
  if(dds_out) {
    check_make_dir("results/")
    saveRDS(dds, paste0("results/", dds_out))
  }
  # Ensure correct format for comparison
  if(typeof(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  # Extract results
  message("Extracting results")
  res <- DESeq2::results(dds, contrast = c(groupCol, comparison))
  
  # Add TPM
  if(tpm & typeof(tpm) == "character"){
    if(!file.exists(tpm)) stop("Database file is missing!\\nLooking for: ", tpm)
    res <- addTPM(res, samples, tpm)
    }
  
  return(res)
}
