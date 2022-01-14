#' Run DESeq2
#' 
#' @import DESeq2
#' @import sva
#' @export
runDESeq2 <- function(txCount,
                      metadata,
                      groupCol,
                      comparison,
                      design,
                      preFilter = 10,
                      parallel = FALSE,
                      cores = 4){
  message("Running DESeq2")
  # Check for group column
  check_colname(colnames(metadata), col_name = groupCol, location = "metadata")
  
  # Register parallel
  if(parallel & is.numeric(cores)) {
    missing_package("BiocParallel", "Bioc")
    BiocParallel::register(BiocParallel::MulticoreParam(4))
  } 
  
  
  ### Ensure rows in metadata matches columns in the count matrix
  txCount <- txCount[, metadata$id]
  
  ### Create DDS from count matrix
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = txCount,
                                        colData = metadata,
                                        design = design) # + condition
  ### Pre-filtering
  if(preFilter){
    message("Pre-filtering with row sum of >=10")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
  }
  
  
  ### Define baseline
  # dds[[groupCol]] <- relevel(dds[[groupCol]], ref = baseline)
  
  ### Run DESeq2
  dds <- DESeq2::DESeq(dds)
  return(dds)
}
