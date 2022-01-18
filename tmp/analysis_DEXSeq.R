# DEXSeq analysis ----

library(DEXSeq)


runDEXseq <- function(dds){
  feat_group <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  svs <- as.character(design(dds))[2] %>% 
    stringr::str_split("\\+ ", n = 2, simplify = TRUE) %>% 
    .[2]
  des <- as.formula(
    paste0("~ sample + exon + condition:exon + ", svs)
  )
  sampleData <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% #row.names = .$id) %>%
    dplyr::select(all_of(groupCol), starts_with("sv")) %>% 
    dplyr::rename(condition = all_of(groupCol))
  #   rename(sample = id)
  # sampleData <- sampleData%>% 
  #   bind_rows(sampleData) %>% 
  #   mutate(sample = factor(sample),
  #          exon = factor(c(rep("this", nrow(sampleData)),
  #                          rep("others", nrow(sampleData))))) %>% 
  #   as.data.frame(row.names = as.character(.$sample))
  
    
  
  # Convert to DEXSeq object
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = counts(dds),
    sampleData = sampleData,
    design = des,
    groupID = feat_group[, 1],
    featureID = feat_group[, 2]
    )
  
}

res_dexseq <- DEXSeq::DEXSeq(dxd, BPPARAM=MulticoreParam(workers=8))

# Shrinkage diagostic
# plotDispEsts( dxd )





