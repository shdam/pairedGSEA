
#' Prepare DEXSeq dataset from DESeq2 dataset
#' 
#' @import DEXSeq
#' @export
prepDEXSeq <- function(dds, groupCol){
  
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
    dplyr::rename(condition = dplyr::all_of(groupCol))
  
  
  ### Convert to DEXSeq object
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = counts(dds),
    sampleData = sampleData,
    design = des,
    groupID = feat_group[, 1],
    featureID = feat_group[, 2]
  )
  
 return(dxd) 
}