#' Run DESeq2 and DEXSeq analyses
#' 
#' @noRd
#' @param row A row from the data.frame of experiments generated with \code{combine_experiments}
run_experiment <- function(row, archs4db = NULL, tx_count = NULL, group_col = "group_nr", tpm = TRUE, prefilter = 10, parallel = TRUE, gtf = NULL, deseq_only = FALSE, store_results = TRUE, run_sva = TRUE){
  
  if(typeof(row) == "character"){ # Convert apply-made row to tibble
    row <- tibble::as_tibble(row, rownames = "names") %>% 
      tidyr::pivot_wider(values_from = value, names_from = names)
  }
  
  
  message("Running on ", row$study)
  
  ### Load metadata
  md_file <- row$filename
  data_name <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove(".csv") 
  
  ### Define tpm file
  if(tpm) tpm <- stringr::str_replace(archs4db, "counts", "tpm")
  ### Define experiment details
  comparison_title <- ifelse(run_sva,
                             yes = row$`comparison_title (empty_if_not_okay)`,
                             no = paste0(row$`comparison_title (empty_if_not_okay)`, "_no_sva"))
  baseline_case <- row$`comparison (baseline_v_condition)` %>% stringr::str_split(pattern = "v", simplify = TRUE) %>% as.character()
  experiment_title <- paste0(data_name, "_", comparison_title)
  
  
  ### Prepare for DE
  if(is.null(tx_count)){
    tx_count <- prepare_tx_count(
      metadata = md_file,
      gtf = gtf,
      archs4db = archs4db,
      group_col = group_col,
      baseline_case = baseline_case
    )
  }
  
  
  results <- pairedGSEA::paired_diff(
    tx_count = tx_count,
    metadata = md_file,
    group_col = group_col,
    sample_col = "id",
    baseline = baseline_case[1],
    case = baseline_case[2],
    experiment_title = experiment_title,
    run_sva = run_sva,
    prefilter = prefilter,
    fit_type = "local",
    store_results = store_results,
    quiet = FALSE,
    deseq_only = deseq_only,
    parallel = parallel,
    BPPARAM = BiocParallel::bpparam()
  )
  
  return(NULL)
  
}


#' Combine experiments
#' @noRd
combine_experiments <- function(md_dir = "metadata"){
  ### List metadata files
  md_files <- list.files(md_dir, full.names = TRUE)
  ### Combine experiments
  experiments <- lapply(md_files, FUN = function(x) {df <- readxl::read_xlsx(x); df$filename <- x; df})  %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(!is.na(`comparison_title (empty_if_not_okay)`))
  return(experiments)
}

#' Load data from ARCHS4 database
#' @noRd
load_archs4 <- function(samples, archs4db, gtf = NULL){
  
  pairedGSEA:::check_missing_package("rhdf5", repo = "Bioc")
  
  # Check that archsdb exists
  if(!file.exists(archs4db) | stringr::str_ends(archs4db, ".h5", negate = TRUE)){
    stop(paste("No ARCHS4 database was found at", archs4db))
  }
  
  ### Extract count data of interest
  # Retrieve information from compressed data
  
  my_ids <- rhdf5::h5read(archs4db, "/meta/samples/geo_accession")
  tx     <- rhdf5::h5read(archs4db, "/meta/transcripts/ensembl_transcript_id")
  if(is.null(gtf)) gene <- rhdf5::h5read(archs4db, "/meta/transcripts/ensembl_gene_id")
  
  # Identify columns to be extracted
  if(!all( samples %in% my_ids )) stop("Some of the chosen samples", samples, "are not in the database.")
  sample_locations <- which(my_ids %in% samples)
  
  # Extract gene expression from compressed data
  tx_count <- archs4db %>% 
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
      dplyr::right_join(dplyr::tibble(transcript = tx), by = "transcript") %>% 
      .[match(tx, .$transcript), ] %>% 
      tidyr::unite(gene, transcript, col = "gene_tx", sep = ":") %>% 
      dplyr::pull("gene_tx")
  }
  
  
  rownames(tx_count) <- gene_tx
  colnames(tx_count) <- my_ids[sample_locations]
  
  return(tx_count)
}

#' Add TPM
#' 
#' @noRd
add_tpm <- function(deseq_results, samples, archs4db_tpm, gtf = NULL){
  ### Add TPM values
  tx_tpm <- load_archs4(samples, archs4db_tpm, gtf = gtf)
  
  deseq_results$tpm <- tx_tpm[rownames(deseq_results), ] %>% 
    rowMeans()
  return(deseq_results)
}





#' Prepare DESeqDataSet from ARCHS4 database
#' @param archs4db The path to an ARCHS4 database file
#' @param gtf Optional: Adds GTF gene and transcript labels to the ARCHS4 data extract
#' @param samples The column in the metadata that specifies the samples
#' @inheritParams paired_gsea
#' @inheritParams prepare_metadata
prepare_tx_count <- function(metadata,
                             group_col,
                             baseline_case,
                             archs4db = NULL,
                             gtf = NULL,
                             samples = "id"){
  # Error tests
  ## Ensure data is provided
  if(is.null(archs4db)) stop("Please provide an Archs4 database.")
  ## Look for database file
  if(!is.null(archs4db)) if(!file.exists(archs4db)) stop("Database file is missing!\nLooking for: ", archs4db)
  
  # Loading metadata
  metadata <- pairedGSEA:::prepare_metadata(metadata, group_col, baseline_case)
  
  # Define samples
  if(samples %in% colnames(metadata)) {samples <- metadata[[samples]]
  } else if(!(typeof(samples) == "character" & length(samples) > 1)) stop("Please specificy 'samples' as a column in metadata or as a vector of samples in database.")
  
  # Load count matrix
  if(!is.null(archs4db)) tx_count <- load_archs4(samples, archs4db, gtf)
  # Ensure rows in metadata matches columns in the count matrix
  tx_count <- tx_count[, samples]
  
  return(tx_count) 
}