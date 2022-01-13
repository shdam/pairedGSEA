# nice /home/ctools/opt/R-4.0.5/bin/R

### Load package
pkgload::load_all()
### Load metadata
md_file <- "tmp/1_GSE154968.xlsx"
metadata <- readxl::read_excel(md_file)


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
# DESeq2 ----

### Define experiment detals
design <- ~ group_nr
baseline <- 2
groupCol <- "group_nr"

### Run DESeq2
dds <- runDESeq2(txCount = txCount,
                 metadata = metadata,
                 groupCol = groupCol,
                 baseline = baseline,
                 design = design,
                 preFilter = 10,
                 parallel = FALSE,
                 cores = 4)
# DESeq2::resultsNames(dds) # lists the coefficients
res <- DESeq2::results(dds, name = resultsNames(dds)[2])

### Add TPM values

txTPM <- loadArchs4(samples, archs4db_tpm)

res$tpm <- txTPM[rownames(res), ] %>% 
  rowSums()


# summary(res)




# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")


# DEXSeq ----

# pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
# list.files(pythonScriptsDir)