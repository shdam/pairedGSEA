# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all()
### Load metadata
md_file <- "tmp/1_GSE154968.xlsx"
metadata <- readxl::read_excel(md_file)

### Define experiment details
comparison <- "2v1"
comparison <- stringr::str_split(comparison, "v", simplify = T)
groupCol <- "group_nr"
metadata[[groupCol]] <- factor(metadata[[groupCol]], levels = comparison)
### Define samples of interest
samples <- metadata$id
### Define file to read from
archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
# TPM file
archs4db_tpm <- archs4db %>% 
  stringr::str_replace("counts", "tpm")
if(!file.exists(archs4db) | !file.exists(archs4db_tpm)) stop("Database file is missing!\\nLooking for: ", archs4db, "and", archs4db_tpm)


### Load count matrix
txCount <- loadArchs4(samples, archs4db)
### Ensure rows in metadata matches columns in the count matrix
txCount <- txCount[, metadata$id]
# Remove low counts
keep <- rowSums(txCount) >= 10
txCount <- txCount[keep,]


### Create DDS from count matrix
dds <- DESeq2::DESeqDataSetFromMatrix(countData = txCount,
                                      colData = metadata,
                                      design = ~1)
# SVA ----
# Normalize counts with DESeq2
normCounts <- DESeq2::normTransform(dds) %>% 
  SummarizedExperiment::assay()

# Define model matrix 
mod1 <- model.matrix(~metadata[[groupCol]])
mod0 <- cbind(mod1[, 1])
# Run SVA
svseq <- sva(normCounts, mod1, mod0)
cat("\\n")
# Store surrogate variables and rename for ease of reference
svs <- tibble::as_tibble(svseq$sv, .name_repair = "minimal")
colnames(svs) <- paste0("sv", 1:svseq$n.sv)
# Add svs to metadata
metadata <- dplyr::bind_cols(metadata, svs)
# Redefine dds colData to metadata
SummarizedExperiment::colData(dds) <- S4Vectors::DataFrame(metadata)
# Redefine design formula to include svs
DESeq2::design(dds) <- as.formula(paste0("~", groupCol, "+", stringr::str_c("sv",1:svseq$n.sv, collapse = "+")))


# DESeq2 ----

### Run DESeq2
dds <- DESeq2::DESeq(dds)

# Extract results
res <- DESeq2::results(dds, contrast = c(groupCol, comparison))


### Add TPM values
txTPM <- loadArchs4(samples, archs4db_tpm)

res$tpm <- txTPM[rownames(res), ] %>% 
  rowSums()

# Summary of results
summary(res)


# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")


# DEXSeq ----

# pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
# list.files(pythonScriptsDir)