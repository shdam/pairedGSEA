# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all()
### Load metadata
md_file <- "tmp/1_GSE154968.xlsx"
### Define file to read from
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
### Define tpm file
tpm <- stringr::str_replace(archs4db, "counts", "tpm")
### Define experiment details
comparison <- "2v1"
groupCol <- "group_nr"

### Run in parallel
BiocParallel::register(BiocParallel::MulticoreParam(4))

### Prepare for DE
dds <- prepDE(md = md_file,
              archs4db = archs4db,
              groupCol = groupCol,
              comparison = comparison,
              prefilter = 10)
### Run DESeq2
res <- runDESeq2(dds,
                 groupCol = groupCol,
                 comparison = comparison,
                 samples = dds$id,
                 tpm = tpm)

# Summary of results
summary(res)


# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

