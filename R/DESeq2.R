#' Run DESeq2
#' 
#' @import DESeq2
#' @export
runDESeq2 <- function(txCount,
                      metadata,
                      group_col,
                      baseline,
                      design,
                      parallel = FALSE,
                      cores = 4){
  # Check for group column
  check_colname(colnames(metadata), col_name = group_col, location = "metadata")
  metadata[[group_col]] <- as.factor(metadata[[group_col]])
  
  # Register parallel
  if(parallel & is.numeric(cores)) {
    missing_package("BiocParallel", "Bioc")
    BiocParallel::register(MulticoreParam(4))
  } 
  
  
  ### Ensure rows in metadata matches columns in the count matrix
  txCount <- txCount[, metadata$id]
  
  ### Create DDS from count matrix
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = txCount,
                                        colData = metadata,
                                        design = design) # + condition
  # Define baseline
  dds[[group_col]] <- relevel(dds[[group_col]], ref = baseline)
  
  ### Run DESeq2
  dds <- DESeq2::DESeq(dds)
  return(dds)
}
