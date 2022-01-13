# DEXSeq analysis ----

library(DEXSeq)

inDir = system.file("extdata", package="pasilla")
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
basename(countFiles)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

# Creating sample table
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3", 
                 "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",  
                "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end", 
               "single-end", "single-end", "paired-end", "paired-end" ) )


dxd <- DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = flattenedFile )

# Normalization
dxd <- estimateSizeFactors( dxd )

# Estimate dispersion
dxd <- estimateDispersions( dxd )

# Shrinkage diagostic
plotDispEsts( dxd )
