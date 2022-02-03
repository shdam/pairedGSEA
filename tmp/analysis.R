# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
### Load metadata
md_file <- "tmp/1_GSE154968.xlsx"
### Define file to read from
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
### Define tpm file
tpm <- stringr::str_replace(archs4db, "counts", "tpm")
### Define experiment details
comparison <- "2v1"
groupCol <- "group_nr"

# Load GTF ----
track <- rtracklayer::import("/home/databases/archs4/v11/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf.gz")
gtf <- tibble::tibble(gene = track$gene_id, transcript = track$transcript_id)
saveRDS(gtf, "gtf.rds")
gtf <- readRDS("gtfextract.rds")

### Run in parallel
# BiocParallel::register(BiocParallel::MulticoreParam(4))

### Prepare for DE ----
dds <- prepDE(md = md_file,
              gtf = gtf,
              archs4db = archs4db,
              groupCol = groupCol,
              comparison = comparison,
              prefilter = 10)

### Prepare DEXSeq
dxd <- prepDEXSeq(dds, groupCol)

### Run DESeq2
res_deseq2 <- runDESeq2(dds,
                        groupCol = groupCol,
                        comparison = comparison,
                        samples = dds$id,
                        tpm = tpm,
                        parallel = TRUE,
                        BPPARAM = BiocParallel::bpparam())#, dds_out = "deseq2_1_GSE154968.RDS")

### Run DEXSeq
res_dexseq <- DEXSeq::DEXSeq(dxd,
                             BPPARAM = BiocParallel::bpparam(),
                             quiet = FALSE)

dataname <- basename(md_file) %>% 
  stringr::str_remove(".xlsx") %>% 
  stringr::str_remove("csv") 

saveRDS(res_deseq2, paste0("results/", dataname, "deseq2res.RDS"))
saveRDS(res_dexseq, paste0("results/", dataname, "dexseqres.RDS"))


plot(metadata(dxr)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(dxr)$lo.fit, col="red")
abline(v=metadata(dxr)$filterTheta)
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")


res <- readRDS("results/32_GSE117523_deseq2res_CDK12 overexpression.RDS")
dxr <- readRDS("results/32_GSE117523_dexseqres_CDK12 overexpression.RDS")

dxr2 <- dxr %>% filter(padj < 0.05)
res2 <- res %>% filter(padj < 0.05)

res2 %>% 
  filter(log2FC > 0) %>% 
  select(gene) %>% 
  rename(Ensembl = gene) %>% 
  readr::write_csv("results/deseq2_genes_up.csv", col_names = TRUE)
res2 %>% 
  filter(log2FC < 0) %>% 
  select(gene) %>% 
  rename(Ensembl = gene) %>% 
  readr::write_csv("results/deseq2_genes_down.csv", col_names = TRUE)
  

# resLFC <- lfcShrink(res=res)

# g:profiler ----

# keep only the significant genes
results_sig = subset(res, padj < 0.05)
# get the significant up-regulated genes
up = subset(results_sig, log2FoldChange > 0)
# get the significant down-regulated genes
down = subset(results_sig, log2FoldChange < 0)
row.names(up) <- row.names(up) %>% str_split(":", simplify = T) %>% .[,2]
row.names(down) <- row.names(down) %>% str_split(":", simplify = T) %>% .[,2]

gp_up = gost(row.names(up), organism = "hsapiens")
gp_down = gost(row.names(down), organism = "hsapiens") 

head(gp_up$result)


up_ordered = up[order(up$log2FoldChange, decreasing = TRUE),]
gp_up_ordered = gost(row.names(up_ordered), organism = "hsapiens",
                     ordered_query = TRUE)
head(gp_up_ordered$result, 8)

gostplot(gp_up, interactive = TRUE)

multi_gp = gost(list("up-regulated" = row.names(up),
                     "down-regulated" = row.names(down)))
p2 = gostplot(multi_gp, interactive = TRUE)



# Convert gene names DESeq2 ----

res_df <- res %>% as_tibble(rownames = "names") %>% 
  tidyr::separate(names, sep = ":", into = c("gene", "transcript"))


results_genes = gconvert(res_df$transcript, organism = "hsapiens", 
                         target = "UNIPROT_GN", filter_na = FALSE)
head(results_genes)

res_df <- res_df %>% 
  arrange(padj)

res_df <- res_df %>% 
  left_join(results_genes, by = c("transcript" = "input"))


# Convert gene names DEXSeq ----

dxr_df <- dxr %>% as_tibble(rownames = "names") %>% 
  tidyr::separate(names, sep = ":", into = c("gene", "transcript"))


dxr_genes = gconvert(dxr_df$transcript, organism = "hsapiens", 
                         target = "UNIPROT_GN", filter_na = FALSE)
head(dxr_genes)

dxr_df <- dxr_df %>% 
  arrange(padj)

dxr_df <- dxr_df %>% 
  left_join(dxr_genes, by = c("transcript" = "input"))


# keep only the significant genes
dxr_sig = subset(dxr, padj < 0.05)
# get the significant up-regulated genes
dxr_up = subset(dxr_sig, log2fold_1_2 > 0)
# get the significant down-regulated genes
dxr_down = subset(dxr_sig, log2fold_1_2 < 0)
row.names(dxr_up) <- row.names(dxr_up) %>% str_split(":", simplify = T) %>% .[,2]
row.names(dxr_down) <- row.names(dxr_down) %>% str_split(":", simplify = T) %>% .[,2]

dxr_gp_up = gost(row.names(dxr_up), organism = "hsapiens")
dxr_gp_down = gost(row.names(dxr_down), organism = "hsapiens") 

head(dxr_gp_up$result)


dxr_up_ordered = dxr_up[order(dxr_up$log2fold_1_2, decreasing = TRUE),]
dxr_gp_up_ordered = gost(row.names(dxr_up_ordered), organism = "hsapiens",
                     ordered_query = TRUE)
head(dxr_gp_up_ordered$result, 8)

gostplot(dxr_gp_up, interactive = TRUE)

dxr_multi_gp = gost(list("up-regulated" = row.names(dxr_up),
                     "down-regulated" = row.names(dxr_down)))
p2 = gostplot(dxr_multi_gp, interactive = TRUE)



# Integrating with external tools for visualisations 

library(clusterProfiler)
library(enrichplot)
library(DOSE) # needed to convert to enrichResult object


up_names = gconvert(row.names(up), target = "UNIPROT_GN")
down_names = gconvert(row.names(down), target = "UNIPROT_GN")

# enrichment analysis using gene names
multi_gp = gost(list("up-regulated" = up_names$name,
                     "down-regulated" = down_names$name), multi_query = FALSE, evcodes = TRUE)

# modify the g:Profiler data frame
gp_mod = multi_gp$result[,c("query", "source", "term_id",
                            "term_name", "p_value", "query_size", 
                            "intersection_size", "term_size", 
                            "effective_domain_size", "intersection")]
gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                  "query_size", "Count", "term_size", "effective_domain_size", 
                  "geneID", "GeneRatio", "BgRatio")
gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
row.names(gp_mod) = stringr::str_c(1:nrow(gp_mod), gp_mod$ID, sep = "_")

# define as compareClusterResult object
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)

# define as enrichResult object
gp_mod_enrich  = new("enrichResult", result = gp_mod)

enrichplot::dotplot(gp_mod_cluster)

barplot(gp_mod_enrich, showCategory = 40, font.size = 16) + 
  ggplot2::facet_grid(~Cluster) +
  ggplot2::ylab("Intersection size")



# AnnotationDb ----

library("AnnotationDbi")
library("org.Hs.eg.db")

res2$symbol <- mapIds(org.Hs.eg.db,
                     keys=res2$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res2$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res) %>% stringr::str_split(":", simplify = T) %>% .[,1],
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
# resOrdered <- res2[order(res2$padj),]
# head(resOrdered)

res %>% 
  as_tibble() %>% 
  filter(log2FC > 0 & padj < 0.05) %>% 
  select(symbol) %>% 
  rename(`Gene symbol` = symbol) %>% 
  readr::write_csv("results/deseq2_symbol_genes_up.csv", col_names = TRUE)

res %>% 
  as_tibble() %>% 
  filter(log2FC < 0 & padj < 0.05) %>% 
  select(symbol) %>% 
  rename(`Gene symbol` = symbol) %>% 
  readr::write_csv("results/deseq2_symbol_genes_down.csv", col_names = TRUE)



dxr2$symbol <- mapIds(org.Hs.eg.db,
                     keys=dxr2$groupID,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

dxrOrdered <- dxr2[order(dxr2$padj),]
resOrdered <- res2[order(res2$padj), ]
# dxrOrdered %>% dplyr::select(symbol, log2FC_baseline_vs_condition, everything())
# resOrdered %>% dplyr::select(symbol, log2FC, everything())

geneID <- dxrOrdered$groupID[1]
dxr %>% as_tibble %>%  filter(groupID == geneID) %>% pull(featureID)
dxrOrdered$symbol[1]

plotDEXSeq2(dxr, geneID, names = TRUE, legend=T, cex.axis=1.2, cex=1.3, lwd=2)



resdxr <- res2 %>%
  full_join(dxrOrdered, by = c("transcript" = "featureID", "symbol"))# %>% 
  # filter(transcript != "missing" & groupID != "missing")

library(ggplot2)
library(tidyverse)
resdxr %>% 
  filter(!is.na(log2FC) & !is.na(log2FC_baseline_vs_condition)) %>% 
  ggplot(aes(x = log2FC, y = log2FC_baseline_vs_condition)) +
  geom_point()
  

resdxr %>% 
  filter(is.na(log2FoldChange) | is.na(log2fold_1_2)) %>% 
  dplyr::select(-gene, -transcript) %>% 
  View

table(table(dxr$groupID))
table(table(res$gene))

dxr2 <- dxr %>% as_tibble(rownames = "transcript") %>% filter(padj < 0.1)

oneGenes <- names(table(dxr$groupID)[which(table(dxr$groupID)==1)])
all(is.na(dxr[which(dxr$groupID %in% oneGenes),"pvalue"]))
res2 %>% 
  filter(gene %in% names(table(dxr$groupID)[which(table(dxr$groupID)==1)]))


plotDEXSeq2( dxr, "missing",
             legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


# Correlation check ----

dxr %>% 
  left_join(res, by= c("featureID" = "transcript")) %>% 
  filter(padj.x < 0.05 & padj.y < 0.05) %>% 
  select(featureID, log2FC, log2FC_baseline_vs_condition)