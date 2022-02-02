
#' Prepare DEXSeq dataset from DESeq2 dataset
#' 
#' @import DEXSeq
#' @export
runDEXSeq <- function(dds, groupCol, comparison){
  message("Initiating DEXSeq")
  # Ensure correct format for comparison
  if(typeof(comparison) != "list") comparison <- stringr::str_split(comparison, "v", simplify = T)
  ### Extract feature and group from rownames of DESeq2 opbject
  feat_group <- rownames(dds) %>% 
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
    dplyr::mutate(condition = dplyr::case_when(condition == comparison[1] ~ 1,
                                               condition == comparison[2] ~ 2) %>% 
                    factor(levels = c(1, 2)))
  
  
  ### Convert to DEXSeq object
  message("Creating DEXSeqDataSet")
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = counts(dds),
    sampleData = sampleData,
    design = des,
    groupID = feat_group[, 1],
    featureID = feat_group[, 2]
  )
  
  ### Run DEXSeq
  message("Running DEXSeq")
  dxr <- DEXSeq::DEXSeq(dxd,
                        reducedModel = as.formula(
                          paste0("~ sample + exon + ", stringr::str_c(svs, ":exon"))
                          ),
                        BPPARAM = BiocParallel::bpparam(),
                        quiet = FALSE)
  dxr <- dxr %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(log2FC_baseline_vs_condition = log2fold_2_1)
  ### Redefine condition to original
  # DEXSeq::sampleAnnotation(dxr)$condition <- dplyr::case_when(DEXSeq::sampleAnnotation(dxr)$condition == "baseline" ~ comparison[1],
  #                                                     TRUE ~ comparison[2]) %>% 
  #   factor(levels = comparison)
  
 return(dxr)
}