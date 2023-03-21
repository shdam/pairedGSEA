#' Paired Over-Representation Analysis
#' 
#' paired_ora uses \code{\link[fgsea:fora]{fora()}} to run the
#' over-representation analysis.
#' First the aggregated pvalues are adjusted using the
#' Benjamini & Hochberg method.
#' The analysis is run on all significant genes found by
#' DESeq2 and DEXSeq individually.
#' I.e., two runs of fora() are executed and subsequently
#' joined into a single object.
#' You can use \code{\link[pairedGSEA:prepare_msigdb]{prepare_msigdb()}}
#' to create a list of gene_sets.
#' 
#' @param paired_diff_result The output of 
#' \code{\link[pairedGSEA:paired_diff]{paired_diff()}}
#' @param gene_sets List of gene sets to analyse
#' @param cutoff (Default: 0.05) Adjusted p-value cutoff for
#' selecting significant genes
#' @param min_size (Default: 25) Minimal size of a gene set to test.
#' All pathways below the threshold are excluded.
#' @param experiment_title Title of your experiment. Your results will be
#' stored in paste0("results/", experiment_title, "_fora.RDS").
#' @param quiet (Default: FALSE) Whether to print messages
#' @param expression_only (Default: FALSE) A logical that indicates whether
#' to only run DESeq2 analysis. Not generally recommended.
#' @importFrom dplyr rename filter arrange rename mutate full_join
#' @importFrom tibble as_tibble
#' @import fgsea
#' @family paired
#' @export 
#' @return A data.frame
#' @usage 
#' paired_ora(
#'     paired_diff_result,
#'     gene_sets,
#'     cutoff = 0.05,
#'     min_size = 25,
#'     experiment_title = NULL,
#'     expression_only = FALSE,
#'     quiet = FALSE
#'     )
#' @examples 
#' data(example_diff_result)
#' gene_sets <- pairedGSEA::prepare_msigdb()
#' 
#' paired_ora(example_diff_result, gene_sets)
#' 
paired_ora <- function(
        paired_diff_result,
        gene_sets,
        cutoff = 0.05,
        min_size = 25,
        experiment_title = NULL,
        expression_only = FALSE,
        quiet = FALSE){
    
    # Check column names are as expected
    if(expression_only) paired_diff_result <- paired_diff_result %>% 
            dplyr::rename(pvalue_expression = pvalue, padj_expression = padj)
    check_colname(paired_diff_result, "pvalue_expression", "paired_diff_result")
    if(!expression_only) check_colname(
        paired_diff_result, "pvalue_splicing", "paired_diff_result")
    check_colname(paired_diff_result, "gene", "paired_diff_result")

    # Significant genes
    if(!quiet)  message("Identifying differentially expressed genes")
    genes_expression <- paired_diff_result %>% 
        dplyr::filter(
            !is.na(pvalue_expression) & !is.na(gene),
            padj_expression < cutoff) %>%
        dplyr::arrange(padj_expression)

    if(!expression_only){
        genes_splicing <- paired_diff_result %>% 
            dplyr::filter(
                !is.na(pvalue_splicing) & !is.na(gene),
                padj_splicing < cutoff) %>%
            dplyr::arrange(padj_splicing)
    }

    # fora
    if(!quiet) message("Running over-representation analysis")
    universe <- unique(paired_diff_result$gene)
    ## ORA on DESeq2 results
    ora_expression <- fgsea::fora(
        gene_sets, genes = genes_expression$gene, 
        universe = universe, minSize = min_size) %>% 
        dplyr::rename(size_geneset = size) %>% 
        dplyr::mutate(
            size_genes = nrow(genes_expression),
            size_universe = length(universe),
            relative_risk =
                (overlap / size_geneset) / (size_genes / size_universe),
            enrichment_score = log2(relative_risk + 0.06))
    if(expression_only){
        if(!is.null(experiment_title)){
            if(!quiet) message("Storing fora results")
            store_result(
                ora_expression, paste0(experiment_title, "_ora.RDS"),
                "ORA on only DESeq2 results", quiet = quiet)
            }
        return(ora_expression)
        }

    ## ORA on DEXSeq results
    ora_splicing <- fgsea::fora(
        gene_sets, genes = genes_splicing$gene,
        universe = universe, minSize = min_size) %>% 
        dplyr::rename(size_geneset = size) %>% 
        dplyr::mutate(
            size_genes = nrow(genes_splicing),
            size_universe = length(universe),
            relative_risk =
                (overlap / size_geneset) / (size_genes / size_universe),
            enrichment_score = log2(relative_risk + 0.06))

    if(!quiet) message("Joining result")
    ora_joined <- ora_expression %>% 
        dplyr::full_join(
            ora_splicing,
            by = c("pathway", "size_geneset", "size_universe"),
            suffix = c("_expression", "_splicing")) %>% 
        dplyr::mutate(
            enrichment_score_shift = 
                relative_risk_splicing - relative_risk_expression,
            enrichment_score_shift = 
                log2(abs(enrichment_score_shift))*sign(enrichment_score_shift),
            experiment = experiment_title) %>% 
        dplyr::arrange(enrichment_score_shift) %>% 
        tibble::as_tibble()

    if(!is.null(experiment_title)){
        if(!quiet) message("Storing fora results")
        store_result(
            ora_joined, paste0(experiment_title, "_ora.RDS"),
            "ORA on both DESeq2 and DEXSeq results", quiet = quiet)
    }
    

    return(ora_joined)
}


#' Load MSigDB and convert to names list of gene sets
#' 
#' This function is wrapper around \code{\link[msigdbr:msigdbr]{msigdbr()}}.
#' Please see their manual for details on its use.
#' The function extracts the gene set name and a user-defined gene id type
#' (Defualt: "ensembl_gene").
#' Please make sure the gene IDs match those from your DE analysis.
#' This function will format the gene sets such that they can be directly
#' used with \code{\link[pairedGSEA:paired_ora]{paired_ora()}}.
#'   
#' @param gene_id_type (Default: "ensemble_gene") The gene ID type to extract.
#' The IDs should match the gene IDs from your DE analysis.
#' @inheritParams msigdbr::msigdbr
#' @note Suggested: importFrom msigdbr msigdbr
#' @return A data.frame
#' @examples 
#' gene_sets <- prepare_msigdb(species = "Homo sapiens")
#' @export
#' @usage 
#' prepare_msigdb(
#'     gene_id_type = "ensembl_gene",
#'     species = "Homo sapiens",
#'     category = "C5",
#'     subcategory = NULL
#'     )
prepare_msigdb <- function(
        gene_id_type = "ensembl_gene",
        species = "Homo sapiens", 
        category = "C5",
        subcategory = NULL
        ){
    check_missing_package("msigdbr")

    gene_sets <- msigdbr::msigdbr(
        species = species,
        category = category,
        subcategory = subcategory
        )
    # Split dataframe based on gene set names
    gene_sets <- gene_sets %>% 
        base::split(x = .[[gene_id_type]], f = .$gs_name)
    return(gene_sets)
}
