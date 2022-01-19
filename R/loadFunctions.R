
#' Load data from ARCHS4 database
#' 
#' @import rhdf5
#' @importFrom stringr str_ends
#' @import magrittr 
#' @importFrom tibble as_tibble tibble
loadArchs4 <- function(samples, archs4db){
  # Check that archsdb exists
  if(!file.exists(archs4db) | stringr::str_ends(archs4db, ".h5", negate = TRUE)){
    stop(paste("No ARCHS4 database was found at", archs4db))
  }
   
  ### Extract count data of interest
  # Retrieve information from compressed data
  myIds <- rhdf5::h5read(archs4db, "/meta/samples/geo_accession")
  gene <- rhdf5::h5read(archs4db, "/meta/transcripts/ensembl_gene_id")
  tx    <- rhdf5::h5read(archs4db, "/meta/transcripts/ensembl_transcript_id")
  gene_tx <- stringr::str_c(gene, tx, sep = ":")
  # Identify columns to be extracted
  if(!all( samples %in% myIds )) stop("Some of the chosen samples", samples, "are not in the database.")
  sample_locations = which(myIds %in% samples)
  
  # Extract gene expression from compressed data
  txCount <- archs4db %>% 
    rhdf5::h5read("data/expression",
           index = list(sample_locations, 1:length(tx))) %>% 
    t()
  # Close file
  rhdf5::H5close()
  rownames(txCount) <- gene_tx
  colnames(txCount) <- myIds[sample_locations]
  
  return(txCount)
  
}

#' Prepare metadata
#' @importFrom stringr str_split
#' @importFrom readxl read_xlsx
#' @importFrom readr read_csv
#' @export
prepMeta <- function(md, groupCol, comparison){
  message("Preparing metadata")
  if(typeof(md[1]) == "character" & length(md) == 1){
    if(stringr::str_ends(md, ".xlsx")) md <- readxl::read_excel(md)
    else if(stringr::str_ends(md, "csv")) md <- readr::read_csv(md)
  }
  
  if(groupCol %!in% colnames(md)) stop("Could not find column", groupCol, "in metadata.")
  ### Ensure comparison is on the right format
  if(typeof(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  ### Remove irrelevant groups
  md <- md[md[[groupCol]] %in% comparison, ]
  ### Add comparison levels to metadata
  md[[groupCol]] <- factor(md[[groupCol]], levels = comparison)
  
  return(md)
}

#' Pre-filter
#' 
#' @export
preFilter <- function(txCount, thres = 10){
  # Remove low counts
  keep <- rowSums(txCount) >= thres
  txCount <- txCount[keep,]
  
  return(txCount)
}


#' Add TPM
#' 
#' @export
addTPM <- function(res, samples, archs4db_tpm){
  ### Add TPM values
  txTPM <- loadArchs4(samples, archs4db_tpm)
  
  res$tpm <- txTPM[rownames(res), ] %>% 
    rowSums()
  return(res)
}