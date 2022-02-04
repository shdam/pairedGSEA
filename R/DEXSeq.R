
#' Prepare DEXSeq dataset from DESeq2 dataset
#' 
#' @import DEXSeq
#' @export
runDEXSeq <- function(dds, groupCol, comparison){
  message("Initiating DEXSeq")
  # Ensure correct format for comparison
  if(length(comparison) == 1) comparison <- stringr::str_split(comparison, "v", simplify = T)
  ### Extract group and feature from rownames of DESeq2 object
  group_feat <- rownames(dds) %>% 
    stringr::str_split(":", simplify = TRUE)
  ### Extract the found surrogate variables
  svs <- as.character(DESeq2::design(dds))[2] %>% 
    stringr::str_split("\\+ ", n = 2, simplify = TRUE) %>% 
    .[2]
  ### Add surrogate variables to DEXSeq design formula
  des <- as.formula(
    paste0("~ sample + exon + condition:exon + ", stringr::str_c(svs, ":exon"))
  )
  
  ### Define sample data based on DESeq2 object
  sampleData <- SummarizedExperiment::colData(dds) %>% 
    as.data.frame(row.names = .$id) %>% #row.names = .$id) %>%
    dplyr::select(dplyr::all_of(groupCol), starts_with("sv")) %>% 
    dplyr::rename(condition = dplyr::all_of(groupCol)) %>% 
    dplyr::mutate(condition = dplyr::case_when(condition == comparison[1] ~ "B",      # baseline 
                                               condition == comparison[2] ~ "C") %>%  # condition
                    factor(levels = c("C", "B")))
  
  
  ### Convert to DEXSeq object
  message("Creating DEXSeqDataSet")
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = counts(dds),
    sampleData = sampleData,
    design = des,
    groupID = group_feat[, 1],
    featureID = group_feat[, 2]
  )
  
  ### Run DEXSeq
  message("Running DEXSeq")
  dxr <- DEXSeq::DEXSeq(dxd,
                        reducedModel = as.formula(
                          paste0("~ sample + exon + ", stringr::str_c(svs, ":exon"))
                          ),
                        BPPARAM = BiocParallel::bpparam(),
                        quiet = FALSE)
  dxr <- dxr %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(log2FC_baseline_vs_condition = log2fold_B_C)
  ### Redefine condition to original
  # DEXSeq::sampleAnnotation(dxr)$condition <- dplyr::case_when(DEXSeq::sampleAnnotation(dxr)$condition == "baseline" ~ comparison[1],
  #                                                     TRUE ~ comparison[2]) %>% 
  #   factor(levels = comparison)
  
 return(dxr)
}





perGeneQValue2 <- function (object,
                            gene = "geneID",
                            p = "pvalue",
                            method = DEXSeq:::perGeneQValueExact) 
{
  # Code adapted from DEXSeq::perGeneQValue
  
  stopifnot(is(object, "DEXSeqResults"))
  wTest <- which(!is.na(object$padj))
  pvals <- object[[p]][wTest]
  geneID <- factor(object[[gene]][wTest])
  geneSplit <- split(seq(along = geneID), geneID)
  pGene <- sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))
  theta <- unique(sort(pGene))
  q <- method(pGene, theta, geneSplit)
  res <- rep(NA_real_, length(pGene))
  res <- q[match(pGene, theta)]
  res <- pmin(1, res)
  names(res) <- names(geneSplit)
  stopifnot(!any(is.na(res)))
  return(res)
}


#' Per gene p value aggregation
#' 
#' @importFrom aggregation lancaster
#' @importFrom purrr when
#' @export
perGenePValue <- function (df,
                           gene = "gene",
                           p = "pvalue",
                           weights = "baseMean",
                           lfc = NULL){
  stopifnot(typeof(df) == "list")
  stopifnot(all(c("padj", gene, p, weights, lfc) %in% colnames(df)))
  
  res <- df %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::rename(pvalue = .data[[p]],
                  ensembl_gene = .data[[gene]]) %>% 
    # Prevent warning from Lancaster
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) %>% 
    dplyr::group_by(ensembl_gene) %>% 
    # p value aggregation
    purrr::when(is.null(lfc) ~ 
                  dplyr::summarise(., pvalue = aggregation::lancaster(pvalue,
                                                                      .data[[weights]])),
                !is.null(lfc) ~ 
                  dplyr::summarise(., pvalue = aggregation::lancaster(pvalue,
                                                                   .data[[weights]]),
                                   lfc = weighted.mean(.data[[lfc]], .data[[weights]]))) %>% 
    dplyr::mutate(pvalue = dplyr::case_when(pvalue < 10e-320 ~ 10e-320,
                                            TRUE ~ pvalue)) 
  return(res)
}
