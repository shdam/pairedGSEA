
#' Load data from ARCHS4 database
#' 
#' @import rhdf5
#' @importFrom stringr str_ends
#' @import magrittr 
#' @importFrom tibble as.tibble
loadArchs4 <- function(samples, archs4db){
  # Check that archsdb exists
  if(!file.exists(archs4db) | stringr::str_ends(archs4db, ".h5", negate = TRUE)){
    stop(paste("No ARCHS4 database was found at", archs4db))
  }
   
  ### Extract count data of interest
  # Retrieve information from compressed data
  myIds <- h5read(archs4db, "meta/samples/geo_accession")
  tx    <- h5read(archs4db, "/meta/transcripts/ensembl_transcript_id")
  
  # Identify columns to be extracted
  message("Extracting:", all( samples %in% myIds ))
  sample_locations = which(myIds %in% samples)
  
  # Extract gene expression from compressed data
  txCount <- archs4db %>% 
    h5read("data/expression",
           index = list(sample_locations, 1:length(tx))) %>% 
    t()
  # Close file
  H5close()
  rownames(txCount) <- tx
  colnames(txCount) <- myIds[sample_locations]
  
  return(as.tibble(txCount))
  
}
