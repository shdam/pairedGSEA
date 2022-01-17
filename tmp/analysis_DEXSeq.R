# DEXSeq analysis ----

library(DEXSeq)



# Creating sample table
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3", 
                 "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",  
                "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end", 
               "single-end", "single-end", "paired-end", "paired-end" ) )

runDEXseq <- function(dds){
  feat_group <- rownames(txCount) %>% 
    stringr::str_split(":", simplify = TRUE)
  # Convert to DEXSeq object
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = assay(dds),
    sampleData = SummarizedExperiment::colData(dds),
    design = design(dds),
    featureID = feat_group[, 1],
    groupID = feat_group[, 2])
}

res_dexseq <- DEXSeq::DEXSeq(dds)

# Normalization
dxd <- DEXSeq::estimateSizeFactors( dxd )

# Estimate dispersion
dxd <- DEXSeq::estimateDispersions( dxd )

res_dexseq <- DEXSeq::DEXSeqResults( dxd )

# Shrinkage diagostic
plotDispEsts( dxd )
