#' Prepare metadata
#' @inheritParams prepDE
prepMeta <- function(md, groupCol, comparison){
  message("Preparing metadata")
  if(typeof(md[1]) == "character" & length(md) == 1){
    if(stringr::str_ends(md, ".xlsx")) md <- readxl::read_excel(md)
    else if(stringr::str_ends(md, ".csv")) md <- readr::read_csv(md)
  } else if("data.frame" %!in% class(md)){stop("Please provide path to a metadata file or a data.frame object.")}
  
  if(groupCol %!in% colnames(md)) stop("Could not find column ", groupCol, " in metadata.")
  ### Ensure comparison is on the right format
  if(typeof(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  ### Remove irrelevant groups
  md <- md[md[[groupCol]] %in% comparison, ]
  ### Add comparison levels to metadata
  md[[groupCol]] <- factor(md[[groupCol]], levels = comparison)
  
  return(md)
}

#' Pre-filter
#' @inheritParams prepDE
preFilter <- function(txCount, thres = 10){
  if(thres < 1) return(txCount)
  # Remove low counts
  keep <- rowSums(txCount) >= thres
  txCount <- txCount[keep,]
  
  return(txCount)
}

#' Run SVA
#' @inheritParams prepDE
runSVA <- function(txCount, metadata, groupCol){
  message("Converting to DESeq object")
  ### Create DDS from count matrix
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = txCount,
                                        colData = metadata,
                                        design = ~1)
  message("Normalizing data")
  # Normalize counts with DESeq2
  normCounts <- DESeq2::normTransform(dds) %>% 
    SummarizedExperiment::assay()
  
  # Define model matrix 
  mod1 <- stats::model.matrix(~metadata[[groupCol]])
  mod0 <- cbind(mod1[, 1])
  # Run SVA
  message("Running SVA")
  svseq <- sva::sva(normCounts, mod1, mod0)
  cat("\n")
  # Store surrogate variables and rename for ease of reference
  svs <- tibble::as_tibble(svseq$sv, .name_repair = "minimal")
  colnames(svs) <- paste0("sv", 1:svseq$n.sv)
  message("Redefining DESeq design formula")
  # Add svs to metadata
  metadata <- dplyr::bind_cols(metadata, svs)
  # Redefine dds colData to metadata
  SummarizedExperiment::colData(dds) <- S4Vectors::DataFrame(metadata)
  # Redefine design formula to include svs
  DESeq2::design(dds) <- as.formula(paste0("~", groupCol, "+", stringr::str_c(colnames(svs), collapse = "+")))
  
  return(dds)
}

#' Prepare DE analysis
#' @param md A metadata file or data frame object
#' @param groupCol The metadata column specifying the what group each sample is associated with
#' @param comparison The comparison to use for this particular experiment. Format example: "1v2"
#' @param archs4db The path to an ARCHS4 database file
#' @param txCount The count matrix of an RNA-seq analysis
#' @param gtf Optional: Adds GTF gene and transcript labels to the ARCHS4 data extract
#' @param samples The column in the metadata that specifies the samples or a vector specifying samples
#' @param prefilter The prefilter value, where rowSums lower than the prefilter value will be removed from the count matrix
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
  # if(!is.null(archs4db)) txCount <- pairedGSEA:::loadArchs4(samples, archs4db, gtf)
  ### Ensure rows in metadata matches columns in the count matrix
  txCount <- txCount[, samples]
  ### Pre-filter
  if(prefilter) txCount <- preFilter(txCount, prefilter)
  
  ### SVA
  dds <- runSVA(txCount, metadata, groupCol)
  
  return(dds) 
}


if(FALSE){
  # Eksempel med den tilsendte metadatafil
  
  md_file <- "72_GSE148800.xlsx"
  comparison <- "1v2"
  groupCol <- "group_nr"
  txCount <- "Your count matrix"
  
  dds <- prepDE(md = md_file,
                groupCol = groupCol,
                comparison = comparison,
                txCount = txCount,
                samples = "id",
                prefilter = 10)
}