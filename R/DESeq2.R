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
                   prefilter = 10){
  
  # Look for database file
  if(!file.exists(archs4db)) stop("Database file is missing!\nLooking for: ", archs4db)
  
  
  # Loading metadata
  metadata <- prepMeta(md, groupCol, comparison)
  
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
#' @import BiocParallel
#' @export
runDESeq2 <- function(dds,
                      groupCol,
                      comparison,
                      fitType = "local",
                      tpm = FALSE,
                      samples = NULL,
                      dds_out = FALSE,
                      parallel = FALSE,
                      BPPARAM = BiocParallel::bpparam()){
  # Register parallel

  message("Running DESeq2")
  dds <- DESeq2::DESeq(dds, parallel = parallel, BPPARAM = BPPARAM,
                       fitType = fitType)
  
  if(typeof(dds_out) == "character") {
    check_make_dir("results/")
    saveRDS(dds, paste0("results/", dds_out))
  }
  # Ensure correct format for comparison
  if(typeof(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  # Extract results
  message("Extracting results")
  res <- DESeq2::results(dds, contrast = c(groupCol, comparison), parallel = parallel, BPPARAM = BPPARAM)
  
  # Add TPM
  if(typeof(tpm) == "character"){
    if(!file.exists(tpm)) stop("Database file is missing!\\nLooking for: ", tpm)
    res <- addTPM(res, samples, tpm)
    }
  
  return(res)
}
