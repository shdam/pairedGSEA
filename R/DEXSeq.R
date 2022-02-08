
#' Prepare DEXSeq dataset from DESeq2 dataset
#' 
#' @import DEXSeq
#' @export
runDEXSeq <- function(dds, groupCol, comparison){
  message("Initiating DEXSeq")
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  ### Extract group and feature from rownames of DESeq2 object
  group_feat <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  ### Extract the found surrogate variables
  svs <- as.character(DESeq2::design(dds))[2] %>% 
    stringr::str_split("\\+ ", n = 2, simplify = TRUE) %>% 
    .[2]
  ### Add surrogate variables to DEXSeq design formula
  des <- as.formula(
    paste0("~ sample + exon + condition:exon + ", stringr::str_c(svs, ":exon"))
  )
  
  ### Define sample data based on DESeq2 object
  sampleData <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% #row.names = .$id) %>%
    dplyr::select(dplyr::all_of(groupCol), starts_with("sv")) %>% 
    dplyr::rename(condition = dplyr::all_of(groupCol)) %>% 
    dplyr::mutate(condition = dplyr::case_when(condition == comparison[1] ~ "B",      # baseline 
                                               condition == comparison[2] ~ "C") %>%  # condition
                    factor(levels = c("C", "B")))
  
  
  ### Convert to DEXSeq object
  message("Creating DEXSeqDataSet")
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = counts(dds),
    sampleData = sampleData,
    design = des,
    groupID = group_feat[, 1],
    featureID = group_feat[, 2]
  )
  
  ### Run DEXSeq
  message("Running DEXSeq")
  dxr <- DEXSeq::DEXSeq(dxd,
                        reducedModel = as.formula(
                          paste0("~ sample + exon + ", stringr::str_c(svs, ":exon"))
                          ),
                        BPPARAM = BiocParallel::bpparam(),
                        quiet = FALSE)
  # Rename LFC column for consistency and human-readability
  dxr <- dxr %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(log2FC_baseline_vs_condition = log2fold_B_C)
  
 return(dxr)
}






