### Load package
pkgload::load_all(path = "/home/projects/shd_pairedGSEA")
# pkgload::load_all()

### List metadata files
md_files <- list.files("metadata", full.names = TRUE)

### Combine each and extract the comparisons to be run
experiments <- lapply(md_files, FUN = function(x) {df <- readxl::read_xlsx(x); df$filename <- x; df})  %>% 
  dplyr::bind_rows() %>% 
  dplyr::filter(!is.na(`comparison_title (empty_if_not_okay)`))


gtf <- readRDS("gtfextract.rds")
gtf <- tibble::tibble(gene = gtf$gene, transcript = gtf$transcript)

BiocParallel::register(MulticoreParam(workers = 10))

# row <- experiments[2,]

### Run experiments ----
runExperiment <- function(row){
  
  row <- tibble::as_tibble(row, rownames = "names") %>% 
    tidyr::pivot_wider(values_from = value, names_from = names)
  
  message("Running on ", row$study)
  
  ### Load metadata
  md_file <- row$filename
  dataname <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove("csv") 
  ### Define file to read from
  archs4db <- "/home/databases/archs4/v11/human_transcript_v11_counts.h5"
  ### Define tpm file
  tpm <- stringr::str_replace(archs4db, "counts", "tpm")
  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)`
  experimentTitle <- row$`comparison_title (empty_if_not_okay)`
  groupCol <- "group_nr"
  
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  
  ### Run in parallel
  # BiocParallel::register(BiocParallel::MulticoreParam(4))
  
  ### Prepare for DE
  dds <- prepDE(md = md_file,
                gtf = gtf,
                archs4db = archs4db,
                groupCol = groupCol,
                comparison = comparison,
                prefilter = 10)
  
  # return(dds)}
  
  ### Run DESeq2
  res <- runDESeq2(dds,
                   groupCol = groupCol,
                   comparison = comparison,
                   samples = dds$id,
                   tpm = tpm,
                   gtf = gtf,
                   parallel = TRUE,
                   fitType = "local",
                   BPPARAM = BiocParallel::bpparam())#, dds_out = "deseq2_1_GSE154968.RDS")
  
  
  saveRDS(res, paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
  rm(res)
  
  ### Run DEXSeq
  dxr <- runDEXSeq(dds, groupCol, comparison)
  
  saveRDS(dxr, paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
  
  message(row$study, " is done.")
  
  
  # or to shrink log fold changes association with condition:
  #res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
  
  
}

apply(experiments[6:nrow(experiments),], 1, runExperiment)


### Analyse experiment results ----

# Load gene names
load('/home/databases/archs4/v11/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.tx.annoation.Rdata')
genes <- ensAnnot %>% dplyr::select(gene_id, gene_name) %>% dplyr::distinct()
rm(ensAnnot)

gene_sets <- msigdbr::msigdbr(category = "C5")

analyseExperiment <- function(row){
  
  row <- tibble::as_tibble(row, rownames = "names") %>% 
    tidyr::pivot_wider(values_from = value, names_from = names)
  
  ### Load metadata
  md_file <- row$filename
  dataname <- basename(md_file) %>% 
    stringr::str_remove(".xlsx") %>% 
    stringr::str_remove("csv") 

  ### Define experiment details
  comparison <- row$`comparison (baseline_v_condition)`
  experimentTitle <- row$`comparison_title (empty_if_not_okay)`
  groupCol <- "group_nr"
  
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  
  
  message("Analysing ", row$study, " ", experimentTitle)
  
  # Load results
  message("Loading results")
  
  res <- readRDS(paste0("results/", dataname, "_deseq2res_", experimentTitle, ".RDS"))
  dxr <- readRDS(paste0("results/", dataname, "_dexseqres_", experimentTitle, ".RDS"))
  
  # p value aggregation
  
  message("Aggregating p values")
  
  dxr_agg <- perGenePValue(dxr, gene = "groupID", weights = "exonBaseMean") %>% 
    dplyr::filter(pvalue < 0.05)
  res_agg <- perGenePValue(res, gene = "gene", weights = "baseMean", lfc = "log2FC") %>% 
    dplyr::filter(pvalue < 0.05)
  
  (comb <- dplyr::full_join(res_agg, dxr_agg, by = "ensembl_gene", suffix = c("_deseq2", "_dexseq")))
  
  # Add gene names
  
  message("Adding gene names")
  
  comb$symbol <- comb %>% 
    dplyr::left_join(genes,
                     by = c("ensembl_gene" = "gene_id"))
  
  message("Storing aggregation result")
  saveRDS(comb, paste0("results/", dataname, "_aggpval_", experimentTitle, ".RDS"))
  
  
  message("Gene set enrichment analysis")
  
  resStats <- comb %>% 
    dplyr::filter(!is.na(pvalue_deseq2)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_deseq2)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene) 
  
  dxrStats <- comb %>% 
    dplyr::filter(!is.na(pvalue_dexseq)) %>% 
    dplyr::mutate(pvalue = -log10(pvalue_dexseq)) %>% 
    dplyr::pull(pvalue, name = ensembl_gene) 
  
  fgseaRes <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = resStats,
                                     BPPARAM = BiocParallel::bpparam()
                                     )
  
  fgseaDxr <- fgsea::fgseaMultilevel(pathways = gene_sets,
                                     stats = dxrStats,
                                     BPPARAM = BiocParallel::bpparam()
  )
  
  
}

apply(experiments, 1, analyseExperiment)

