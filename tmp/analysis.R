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

metadata <- prepMeta(md_file, groupCol, comparison)

### Define samples of interest
samples <- metadata$id

# TPM file
archs4db_tpm <- archs4db %>% 
  stringr::str_replace("counts", "tpm")
if(!file.exists(archs4db) | !file.exists(archs4db_tpm)) stop("Database file is missing!\\nLooking for: ", archs4db, "and", archs4db_tpm)


### Load count matrix
txCount <- loadArchs4(samples, archs4db)
### Ensure rows in metadata matches columns in the count matrix
txCount <- txCount[, metadata$id] %>% 
  preFilter()


# SVA ----

dds <- runSVA(txCount, metadata, groupCol)

# DESeq2 ----
res <- runDESeq2(dds, groupCol, comparison, parallel = TRUE)

### Add TPM to results
res <- addTPM(res, samples, archs4db_tpm)

# Summary of results
summary(res)


# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

