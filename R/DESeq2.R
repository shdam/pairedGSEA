#' Pepare DE analysis
#' 
#' @import DESeq2
#' @import sva
#' @export
prepDE <- function(md,
                   groupCol,
                   comparison,
                   archs4db = NULL,
                   txCount = NULL,
                   gtf = NULL,
                   samples = "id",
                   prefilter = 10){
  ### Error tests
  # Esnure data is provided
  if(is.null(archs4db) & is.null(txCount)) stop("Please provide a transcript count matrix or a Archs4 database.")
  # Look for database file
  if(!is.null(archs4db)) if(!file.exists(archs4db)) stop("Database file is missing!\nLooking for: ", archs4db)
  
  # Loading metadata
  metadata <- prepMeta(md, groupCol, comparison)
  
  # Define samples
  if(samples %in% colnames(metadata)) {samples <- metadata[[samples]]
  } else if(!(typeof(samples) == "character" & length(samples) > 1)) stop("Please specificy 'samples' as a column in metadata or as a vector of samples in database.")
  
  ### Load count matrix
  if(!is.null(archs4db)) txCount <- loadArchs4(samples, archs4db, gtf)
  ### Ensure rows in metadata matches columns in the count matrix
  txCount <- txCount[, samples]
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
                      gtf = NULL,
                      dds_out = FALSE,
                      parallel = FALSE,
                      BPPARAM = BiocParallel::bpparam()){
  # Register parallel

  message("Running DESeq2")
  dds <- DESeq2::DESeq(dds, parallel = parallel, BPPARAM = BPPARAM,
                       fitType = fitType, quiet = FALSE)
  
  if(typeof(dds_out) == "character") {
    check_make_dir("results/")
    saveRDS(dds, paste0("results/", dds_out))
  }
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  # Extract results
  message("Extracting results")
  res <- DESeq2::results(dds, contrast = c(groupCol, comparison), parallel = parallel, BPPARAM = BPPARAM)
  
  # Add TPM
  if(typeof(tpm) == "character"){
    if(!file.exists(tpm)) stop("Database file is missing!\\nLooking for: ", tpm)
    res <- addTPM(res, samples, tpm, gtf)
    }
  
  # Convert result to tibble
  res <- res %>% 
    tibble::as_tibble(rownames = "gene_tx") %>% 
    tidyr::separate(gene_tx, into = c("gene", "transcript"), sep = ":") %>% 
    dplyr::rename(log2FC = log2FoldChange)
  
  return(res)
}
