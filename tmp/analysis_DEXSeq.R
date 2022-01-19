# DEXSeq analysis ----

library(DEXSeq)


dxd <- prepDEXSeq(dds, "group_nr")

res_dexseq <- DEXSeq::DEXSeq(dxd, BPPARAM=bpparam())

# Shrinkage diagostic
# plotDispEsts( dxd )





