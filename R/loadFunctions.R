
#' Combine experiments
#' 
#' @importFrom readxl read_xlsx
#' @export
combineExperiments <- function(md_dir = "metadata"){
  ### List metadata files
  md_files <- list.files(md_dir, full.names = TRUE)
  ### Combine experiments
  experiments <- lapply(md_files, FUN = function(x) {df <- readxl::read_xlsx(x); df$filename <- x; df})  %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(!is.na(`comparison_title (empty_if_not_okay)`))
  return(experiments)
}



#' Load data from ARCHS4 database
#' 
#' @import rhdf5
#' @importFrom stringr str_ends
#' @import magrittr 
#' @importFrom tibble as_tibble tibble
loadArchs4 <- function(samples, archs4db, gtf = NULL){
  # Check that archsdb exists
  if(!file.exists(archs4db) | stringr::str_ends(archs4db, ".h5", negate = TRUE)){
    stop(paste("No ARCHS4 database was found at", archs4db))
  }
   
  ### Extract count data of interest
  # Retrieve information from compressed data
  
  myIds <- rhdf5::h5read(archs4db, "/meta/samples/geo_accession")
  tx    <- rhdf5::h5read(archs4db, "/meta/transcripts/ensembl_transcript_id")
  if(is.null(gtf)) gene <- rhdf5::h5read(archs4db, "/meta/transcripts/ensembl_gene_id")
  
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
  
  if(is.null(gtf)){
    gene_tx <- stringr::str_c(gene, tx, sep = ":")
  } else{
    gene_tx <- gtf %>% 
      dplyr::distinct() %>% 
      dplyr::right_join(tibble::tibble(transcript = tx), by = "transcript") %>% 
      .[match(tx, .$transcript), ] %>% 
      tidyr::unite(gene, transcript, col = "gene_tx", sep = ":") %>% 
      dplyr::pull("gene_tx")
  }
  
  
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
  if(thres < 1) return(txCount)
  # Remove low counts
  keep <- rowSums(txCount) >= thres
  txCount <- txCount[keep,]
  
  return(txCount)
}


#' Add TPM
#' 
#' @export
addTPM <- function(res, samples, archs4db_tpm, gtf = NULL){
  ### Add TPM values
  txTPM <- loadArchs4(samples, archs4db_tpm, gtf = gtf)
  
  res$tpm <- txTPM[rownames(res), ] %>% 
    rowMeans()
  return(res)
}



loadGTF <- function(samples, gtf){
  
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



#' Load MSigDB and convert to names list of gene sets
#' 
#' @import msigdbr
#' @importFrom purrr map
#' @export
prepMsigdb <- function(category = "C5"){
  gene_sets <- msigdbr::msigdbr(category = "C5")
  
  gene_sets <- gene_sets %>% 
    dplyr::select(gs_name, ensembl_gene) %>% 
    # Split dataframe based on gene set names
    base::split(.$gs_name) %>% 
    # In each list, extract only the ensemble IDs
    purrr::map(.f = ~ .x$ensembl_gene)
  
  return(gene_sets)
}