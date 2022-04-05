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