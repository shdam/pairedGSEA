library(rhdf5)

### Define samples of interest
samples <- c("GSM1383738","GSM1383739","GSM1383740","GSM1383741","GSM1383742")
### Define file to read from
archsdb <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"

#' Load data from ARCHS4 database
#' 
#' @import rhdf5
#' @importFrom stringr str_
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
