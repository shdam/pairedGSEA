#' Run SVA on DESeqDataSet
#' @inheritParams paired_gsea
run_sva <- function(dds, group_col, quiet = FALSE){
  
  
  if(!quiet) message("Running SVA")
  # Normalize counts with DESeq2 for SVA
  normalized_counts <- DESeq2::normTransform(dds) %>% 
    SummarizedExperiment::assay()

  # Extract metadata
  metadata <- SummarizedExperiment::colData(dds)
  # Define model matrix 
  mod1 <- stats::model.matrix(~metadata[[group_col]])
  mod0 <- cbind(mod1[, 1])
  
  # Run SVA
  svseq <- sva::sva(normalized_counts, mod1, mod0)
  message("Found ", svseq$n.sv, " surrogate variables")
  # Store surrogate variables and rename for ease of reference
  svs <- tibble::as_tibble(svseq$sv, .name_repair = "minimal")
  colnames(svs) <- paste0("sv", 1:svseq$n.sv)
  
  if(!quiet) message("Redefining DESeq design formula\n")
  # Add svs to dds colData
  SummarizedExperiment::colData(dds) <- cbind(metadata, svs)
  # Redefine design formula to include svs
  DESeq2::design(dds) <- as.formula(paste0("~", group_col, "+", stringr::str_c(colnames(svs), collapse = "+")))

  
  return(dds)
}