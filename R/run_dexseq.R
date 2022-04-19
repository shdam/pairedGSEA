#' Prepare DEXSeqDataSet from DESeqDataSet
#' 
#' @inheritParams paired_gsea
run_dexseq <- function(dds,
                       group_col,
                       comparison,
                       dxd_out = FALSE,
                       quiet = FALSE,
                       parallel = FALSE,
                       BPPARAM = BiocParallel::bpparam()){
  
  if(parallel) pairedGSEA:::check_missing_package("BiocParallel", "Bioc")
  
  if(!quiet) message("Initiating DEXSeq")
  # Ensure correct format for comparison
  comparison <- pairedGSEA:::check_comparison(comparison)
  
  # Extract group and feature from rownames of DESeq2 object
  group_feat <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  
  # Extract the found surrogate variables
  svs <- as.character(DESeq2::design(dds))[2] %>% 
    stringr::str_split("\\+ ", n = 2, simplify = TRUE) %>% 
    .[2]
  
  # Add surrogate variables to DEXSeq design formula
  design_formula <- as.formula(
    paste0("~ sample + exon + condition:exon + ", stringr::str_c(svs, ":exon"))
  )
  
  # Define sample data based on DESeq2 object
  sample_data <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% 
    dplyr::select(dplyr::all_of(group_col), dplyr::starts_with("sv")) %>% 
    dplyr::rename(condition = dplyr::all_of(group_col)) %>% 
    dplyr::mutate(condition = dplyr::case_when(condition == comparison[1] ~ "B",      # baseline 
                                               condition == comparison[2] ~ "C") %>%  # condition
                    factor(levels = c("C", "B")))
  
  
  # Convert to DEXSeqDataSet
  if(!quiet) message("Creating DEXSeqDataSet")
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = DESeq2::counts(dds),
    sampleData = sample_data,
    design = design_formula,
    groupID = group_feat[, 1],
    featureID = group_feat[, 2]
  )
  
  # Store DEXSeqDataSet with DEXSeq analysis
  if(typeof(dxd_out) == "character") {
    pairedGSEA:::store_result(dxd, dxd_out, "DEXSeqDataSet", quiet = quiet)
  }
  ### Run DEXSeq
  if(!quiet) message("Running DEXSeq")
  dxr <- DEXSeq::DEXSeq(dxd,
                        reducedModel = as.formula(
                          paste0("~ sample + exon + ", stringr::str_c(svs, ":exon"))
                          ),
                        BPPARAM = ifelse(parallel, BiocParallel::bpparam(), NULL),
                        quiet = quiet)
  # Rename LFC column for consistency and human-readability
  dxr <- dxr %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(log2FC_dexseq = log2fold_B_C)
  
 return(dxr)
}






