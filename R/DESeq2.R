#' Run DESeq2
#' 
#' @import DESeq2
#' @import sva
#' @export
runDESeq2 <- function(dds,
                      groupCol,
                      comparison,
                      parallel = FALSE,
                      cores = 4){
  message("Running DESeq2")
  
  # Register parallel
  if(parallel & is.numeric(cores)) {
    missing_package("BiocParallel", "Bioc")
    BiocParallel::register(BiocParallel::MulticoreParam(4))
  } 
  
  
  ## Ensure comparison is on the right format
  if(type(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  ### Run DESeq2
  dds <- DESeq2::DESeq(dds, parallel = parallel, BPPARAM = cores)
  
  # Extract results
  res <- DESeq2::results(dds, contrast = c(groupCol, comparison))
  
  return(res)
}
