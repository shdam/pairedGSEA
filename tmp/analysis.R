# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all()
### Load metadata
md_file <- "tmp/1_GSE154968.xlsx"
metadata <- readxl::read_excel(md_file)

### Define experiment detals
comparison <- "2v1"
groupCol <- "group_nr"
metadata[[groupCol]] <- factor(metadata[[groupCol]], levels = stringr::str_split(comparison, "v", simplify = T))
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
# Remove low counts
keep <- rowSums(txCount) >= 10
txCount <- txCount[keep,]

# SVA ----

mod1 <- model.matrix(~metadata[[groupCol]])
mod0 <- cbind(mod1[, 1])
svseq <- svaseq(txCount, mod1, mod0)

mod1sv <- cbind(mod1, svseq$sv)
mod0sv <- cbind(mod0, svseq$sv)
# pValuesSv <- f.pvalue(txCount, mod1sv, mod0sv)
# qValuesSv <- p.adjust(pValuesSv, method = "BH")


# DESeq2 ----

colnames(mod1sv) <- c("Intercept", groupCol, 1:svseq$n.sv)
mod1sv
design <- mod1sv

### Run DESeq2
dds <- runDESeq2(txCount = txCount,
                 metadata = metadata,
                 groupCol = groupCol,
                 comparison = comparison,
                 design = design,
                 preFilter = FALSE,
                 parallel = TRUE,
                 cores = 4)
# DESeq2::resultsNames(dds) # lists the coefficients
# contrast <- c(groupCol,
#               rev(stringr::str_split(comparison, "v", simplify = T)))
contrast <- list("add" = resultsNames(dds)[2], "subtract" = resultsNames(dds)[3:(2+svseq$n.sv)])
print("contrast:", contrast)
res <- DESeq2::results(dds, contrast = contrast)


### Add TPM values
txTPM <- loadArchs4(samples, archs4db_tpm)

res$tpm <- txTPM[rownames(res), ] %>% 
  rowSums()


summary(res)


# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")


# DEXSeq ----

# pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
# list.files(pythonScriptsDir)