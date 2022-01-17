# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all()
### Load metadata
md_file <- "tmp/1_GSE154968.xlsx"

### Define file to read from
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"

### Define experiment details
comparison <- "2v1"
# comparison <- stringr::str_split(comparison, "v", simplify = T)
groupCol <- "group_nr"
# metadata[[groupCol]] <- factor(metadata[[groupCol]], levels = comparison)

res <- runDESeq2(md = md_file,
                 archs4db =  archs4db,
                 groupCol = groupCol,
                 comparison = comparison,
                 prefilter = 10,
                 samples = "id",
                 tpm = TRUE,
                 parallel = TRUE)


# Summary of results
summary(res)


# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

