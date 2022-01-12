# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all()

### Load metadata
md <- readxl::read_excel("tmp/1_GSE154968.xlsx")

### Define samples of interest
samples <- emeta$id
### Define file to read from
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"

### Load count matrix
txCount <- loadArchs4(samples, archs4db)


dds <- DESeqDataSetFromMatrix(countData = txCount,
                              colData = md,
                              design= ~ group_nr) # + condition
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
