


res <- readRDS("results/109_GSE171358_deseq2res_Cabozantinib treatment.RDS")
dxr <- readRDS("results/109_GSE171358_dexseqres_Cabozantinib treatment.RDS")
res

qplot(res$pvalue, main = "deseq2 pval")
qplot(dxr$pvalue, main = "dexseq pval")

sum(res$padj < 0.05, na.rm = TRUE)
sum(dxr$padj < 0.05, na.rm = TRUE)

comb <- readRDS("results/109_GSE171358_aggpval_Cabozantinib treatment.RDS")


sum(comb$pvalue_deseq2 < 0.05, na.rm = T)
sum(comb$pvalue_dexseq < 0.05, na.rm = T)

qplot(comb$pvalue_dexseq, main = "dexseq pval")
qplot(comb$pvalue_deseq2, main = "deseq2 pval")

dxr_agg <- perGenePValue(dxr, gene = "groupID", weights = "exonBaseMean", lfc = "log2FC_baseline_vs_condition", type = "dexseq")# %>% 
# filter(pvalue < 0.1)
res_agg <- perGenePValue(res, gene = "gene", weights = "baseMean", lfc = "log2FC", type = "deseq2")# %>% 
# filter(pvalue < 0.1)

(comb <- dplyr::full_join(res_agg, dxr_agg, by = "ensembl_gene", suffix = c("_deseq2", "_dexseq")))


resGenes <- comb %>% 
  filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>%
  mutate(padj = p.adjust(pvalue_deseq2, "fdr")) %>% 
  filter(padj < 0.05) %>% 
  arrange(padj)
dxrGenes <- comb %>% 
  filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>%
  mutate(padj = p.adjust(pvalue_dexseq, "fdr")) %>% 
  filter(padj < 0.05) %>% 
  arrange(padj)

nrow(resGenes)
nrow(dxrGenes)


forares <- fgsea::fora(gene_sets, genes = resGenes$ensembl_gene, 
                       universe = unique(comb$ensembl_gene), minSize = 25)
foradxr <- fgsea::fora(gene_sets, genes = dxrGenes$ensembl_gene,
                       universe = unique(comb$ensembl_gene), minSize = 25)
foraresdxr <- fgsea::fora(gene_sets, genes = intersect(dxrGenes$ensembl_gene, resGenes$ensembl_gene,),
                       universe = unique(comb$ensembl_gene), minSize = 25)

saveRDS(forares, paste0("results/", dataname, "_forares_", experimentTitle, ".RDS"))
saveRDS(foradxr, paste0("results/", dataname, "_foradxr_", experimentTitle, ".RDS"))
saveRDS(foraresdxr, paste0("results/", dataname, "_foraresdxr_", experimentTitle, ".RDS"))

sum(forares$padj < 0.05)
sum(foradxr$padj < 0.05)


qplot(forares$pval)
qplot(foradxr$pval)

resStats <- comb %>% 
  dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = -log10(pvalue_deseq2) * sign(lfc_deseq2)) %>% 
  dplyr::pull(pvalue, name = ensembl_gene)
resStats2 <- comb %>% 
  dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = -log10(pvalue_deseq2)) %>% 
  dplyr::pull(pvalue, name = ensembl_gene)

resStatss <- scale( rank( resStats2))
resStatss2 <- comb %>% 
  dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = stat2) %>% 
  dplyr::pull(pvalue, name = ensembl_gene)

resStatss <- comb %>% 
  dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = stat_deseq2) %>% 
  dplyr::pull(pvalue, name = ensembl_gene)
resStatss2 <- comb %>% 
  dplyr::filter(!is.na(pvalue_deseq2) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = stat2) %>% 
  dplyr::pull(pvalue, name = ensembl_gene)

fgseaRes_oristd <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                          stats = resStats,
                                          scoreType = "std",
                                          eps = 10e-320
)
fgseaRes_oripos <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                          stats = resStats,
                                          scoreType = "pos",
                                          eps = 10e-320
)
fgseaRes_s <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                          stats = resStatss,
                                          scoreType = "std",
                                          eps = 10e-320
)
fgseaRes_s2 <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                          stats = resStatss2,
                                          scoreType = "std",
                                          eps = 10e-320
)

fgseaDxr <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = dxrStats,
                                     scoreType = "pos",
                                   minSize = 25,
                                     eps = 10e-320
)
fgsea <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                      stats = resStatss2,
                                      scoreType = "std",
                                      eps = 10e-320
)

dxrStats <- comb %>% 
  dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = stat_dexseq) %>% 
  dplyr::pull(pvalue, name = ensembl_gene) 
dxrStats2 <- comb %>% 
  dplyr::filter(!is.na(pvalue_dexseq) & !is.na(ensembl_gene)) %>% 
  dplyr::mutate(pvalue = stat_dexseq * sign(lfc_dexseq)) %>% 
  dplyr::pull(pvalue, name = ensembl_gene) 

qplot(resStats, bins = 100)
qplot(dxrStats)
qplot(dxrStats2)

fgseaDxr_s <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                      stats = dxrStats,
                                      scoreType = "pos",
                                      eps = 10e-320
)
fgseaDxr_s2 <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = dxrStats2,
                                     scoreType = "std",
                                     eps = 10e-320
)

fgseaRes <- readRDS("results/109_GSE171358_fgseaRes_Cabozantinib treatment.RDS")
fgseaDxr <- readRDS("results/109_GSE171358_fgseaDxr_Cabozantinib treatment.RDS")

qplot(fgseaRes$pval)
qplot(fgseaDxr$pval)
qplot(fgseaDxr2$pval)
qplot(fgseaDxr3$pval)
qplot(fgseaRes$padj)
qplot(fgseaDxr$padj)


length(intersect(fgseaRes_pos$pathway[which(fgseaRes_pos$padj < 0.05)], fgseaRes_std$pathway[which(fgseaRes_std$padj < 0.05)]))
length(intersect(fgseaDxr_pos$pathway[which(fgseaDxr_pos$padj < 0.05)], fgseaDxr_std$pathway[which(fgseaDxr_std$padj < 0.05)]))

qplot(-log10(dxr$pvalue), -log10(res$pvalue))
qplot(-log10(fgseaDxr2$pval), -log10(fgseaRes$pval))

sum(fgseaRes$padj < 0.05, na.rm = T)
# sum(fgseaRes3$padj < 0.05, na.rm = T)
sum(fgseaDxr$padj < 0.05, na.rm = T)
sum(fgseaDxr2$padj < 0.05, na.rm = T)

dxrtrans <- fgseaDxr2$pathway[which(fgseaDxr2$padj < 0.05)]
restrans <- fgseaRes$pathway[which(fgseaRes$padj < 0.05)]

restrans <- res$transcript %in% res$transcript[res$padj < 0.05]
dxrtrans <- dxr$featureID %in% dxr$featureID[dxr$padj < 0.05]

restrans <- res$transcript[which(res$padj < 0.05)]
dxrtrans <- dxr$featureID[which(dxr$padj < 0.05)]

proxy::simil(restrans, dxrtrans, by_rows = FALSE, method = "Simpson")


comb %>% 
  filter(!is.na(pvalue_deseq2) & !is.na(pvalue_dexseq)) %>% 
  ggplot(aes(x = pvalue_deseq2, y = pvalue_dexseq)) +
  geom_point()
