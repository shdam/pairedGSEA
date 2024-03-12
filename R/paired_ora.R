#' Paired Over-Representation Analysis
#' 
#' paired_ora uses \code{\link[fgsea:fora]{fora}} to run the
#' over-representation analysis.
#' First the aggregated pvalues are adjusted using the
#' Benjamini & Hochberg method.
#' The analysis is run on all significant genes found by
#' \code{\link[DESeq2:DESeq]{DESeq2}} and 
#' \code{\link[DEXSeq:DEXSeq]{DEXSeq}} individually.
#' I.e., two runs of \code{\link[fgsea:fora]{fora}}
#' are executed and subsequently
#' joined into a single object.
#' You can use \code{\link[pairedGSEA:prepare_msigdb]{prepare_msigdb}}
#' to create a list of gene_sets.
#' 
#' @param paired_diff_result The output of 
#' \code{\link[pairedGSEA:paired_diff]{paired_diff}}
#' @param gene_sets List of gene sets to analyse
#' @param cutoff (Default: \code{0.05}) Adjusted p-value cutoff for
#' selecting significant genes
#' @param min_size (Default: \code{25}) Minimal size of a gene set to test.
#' All pathways below the threshold are excluded.
#' @param experiment_title Title of your experiment. Your results will be
#' stored in \code{paste0("results/", experiment_title, "_fora.RDS")}.
#' @param quiet (Default: \code{FALSE}) Whether to print messages
#' @param expression_only (Default: \code{FALSE})
#' A logical that indicates whether
#' to only run \code{\link[DESeq2:DESeq]{DESeq2}} analysis.
#' Not generally recommended.
#' @importFrom S4Vectors DataFrame complete.cases
#' @import fgsea
#' @family paired
#' @export 
#' @return A data.table of merged ORA results
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
#' data("example_diff_result")
#' data("example_gene_sets")
#' 
#' ora <- paired_ora(
#'     example_diff_result,
#'     example_gene_sets)
#' 
#' 
paired_ora <- function(
        paired_diff_result,
        gene_sets,
        cutoff = 0.05,
        min_size = 25,
        experiment_title = NULL,
        expression_only = FALSE,
        quiet = FALSE){
    ## Initial error checks
    stopifnot(is(c(quiet, expression_only), "logical"))
    stopifnot(is(c(cutoff, min_size), "numeric"))
    stopifnot(
        is(paired_diff_result, "data.frame") | is(paired_diff_result, "DFrame")
        )
    
    # Check column names are as expected
    if(expression_only) {
        colnames(paired_diff_result)[
            which(colnames(paired_diff_result) == "pvalue")] <- 
            "pvalue_expression"
        colnames(paired_diff_result)[
            which(colnames(paired_diff_result) == "padj")] <-
            "padj_expression"
    } else{
        check_colname(
            paired_diff_result, "pvalue_splicing", "paired_diff_result")
    }
    check_colname(
        paired_diff_result, "pvalue_expression", "paired_diff_result")
    check_colname(paired_diff_result, "gene", "paired_diff_result")
    
    if(!quiet) message("Running over-representation analyses")
    # fora on expression
    ora_expression <- run_ora(
        paired_diff_result, gene_sets = gene_sets,
        type = "expression", cutoff = cutoff, min_size = min_size
    )
    
    if(expression_only){
        if(!is.null(experiment_title)){
            if(!quiet) message("Storing fora results")
            store_result(
                ora_expression,
                paste0(experiment_title, "_ora.RDS"),
                "ORA on only DESeq2 results", quiet = quiet)
        }
        return(ora_expression)
    }
    # fora on splicing
    ora_splicing <- run_ora(
        paired_diff_result, gene_sets = gene_sets,
        type = "splicing", cutoff = cutoff, min_size = min_size
    )
    
    # fora on paired
    ora_paired <- run_ora(
        paired_diff_result, gene_sets = gene_sets,
        type = "paired", cutoff = cutoff, min_size = min_size
    )
    
    if(!quiet) message("Joining result")
    ora_joined <- join_results(ora_expression, ora_splicing, ora_paired)
    
    if(!is.null(experiment_title)){
        if(!quiet) message("Storing fora results")
        ora_joined$experiment <- experiment_title
        store_result(
            ora_joined, paste0(experiment_title, "_ora.RDS"),
            "ORA on differential analysis results", quiet = quiet)
    }
    return(S4Vectors::DataFrame(ora_joined))
}

#' Paired Functional Class Scoring
#' 
#' paired_fcs uses \code{\link[fgsea:fgseaMultilevel]{fgseaMultilevel}} to run the
#' FCS analysis.
#' First the aggregated pvalues are adjusted using the
#' Benjamini & Hochberg method.
#' The analysis is run on all significant genes found by
#' \code{\link[DESeq2:DESeq]{DESeq2}} and 
#' \code{\link[DEXSeq:DEXSeq]{DEXSeq}} individually.
#' I.e., two runs of \code{\link[fgsea:fgseaMultilevel]{fgseaMultilevel}}
#' are executed and subsequently
#' joined into a single object.
#' You can use \code{\link[pairedGSEA:prepare_msigdb]{prepare_msigdb}}
#' to create a list of gene_sets.
#' 
#' @inheritParams paired_ora
#' @importFrom S4Vectors DataFrame complete.cases
#' @import fgsea
#' @family paired
#' @export 
#' @return A data.table of merged FCS results
#' @examples 
#' data("example_diff_result")
#' data("example_gene_sets")
#' 
#' ora <- paired_fcs(
#'     example_diff_result,
#'     example_gene_sets)
#' 
#' 
paired_fcs <- function(
        paired_diff_result,
        gene_sets,
        min_size = 25,
        experiment_title = NULL,
        expression_only = FALSE,
        quiet = FALSE){
    ## Initial error checks
    stopifnot(is(c(quiet, expression_only), "logical"))
    stopifnot(is(min_size, "numeric"))
    stopifnot(
        is(paired_diff_result, "data.frame") | is(paired_diff_result, "DFrame")
    )
    
    # Check column names are as expected
    if(expression_only) {
        colnames(paired_diff_result)[
            which(colnames(paired_diff_result) == "pvalue")] <- 
            "pvalue_expression"
        colnames(paired_diff_result)[
            which(colnames(paired_diff_result) == "padj")] <-
            "padj_expression"
    } else{
        check_colname(
            paired_diff_result, "pvalue_splicing", "paired_diff_result")
    }
    check_colname(
        paired_diff_result, "pvalue_expression", "paired_diff_result")
    check_colname(paired_diff_result, "gene", "paired_diff_result")
    
    if(!quiet) message("Running functional class scoring analyses")
    # fora on expression
    expression <- run_fcs(
        paired_diff_result, gene_sets = gene_sets,
        type = "expression", min_size = min_size
    )
    
    if(expression_only){
        if(!is.null(experiment_title)){
            if(!quiet) message("Storing fora results")
            store_result(
                ora_expression,
                paste0(experiment_title, "_fcs.RDS"),
                "FCS on only DESeq2 results", quiet = quiet)
        }
        return(expression)
    }
    # fcs on splicing
    splicing <- run_fcs(
        paired_diff_result, gene_sets = gene_sets,
        type = "splicing", min_size = min_size
    )
    
    # fcs on paired
    paired <- run_fcs(
        paired_diff_result, gene_sets = gene_sets,
        type = "paired", min_size = min_size
    )
    
    if(!quiet) message("Joining result")
    joined <- join_results(expression, splicing, paired, type = "fcs")
    
    if(!is.null(experiment_title)){
        if(!quiet) message("Storing fcs results")
        joined$experiment <- experiment_title
        store_result(
            joined, paste0(experiment_title, "_fcs.RDS"),
            "FCS on differential analysis results", quiet = quiet)
    }
    return(S4Vectors::DataFrame(joined))
}

#' Load MSigDB and convert to names list of gene sets
#' 
#' This function is wrapper around \code{\link[msigdbr:msigdbr]{msigdbr()}}.
#' Please see their manual for details on its use.
#' The function extracts the gene set name and a user-defined gene id type
#' (Default: "ensembl_gene").
#' Please make sure the gene IDs match those from your DE analysis.
#' This function will format the gene sets such that they can be directly
#' used with \code{\link[pairedGSEA:paired_ora]{paired_ora()}}.
#'   
#' @param gene_id_type (Default: "ensemble_gene") The gene ID type to extract.
#' The IDs should match the gene IDs from your DE analysis.
#' @inheritParams msigdbr::msigdbr
#' @note Suggested: importFrom msigdbr msigdbr
#' @return A list of gene sets
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
    gene_sets <- base::split(
        gene_sets,
        x = gene_sets[[gene_id_type]],
        f = gene_sets$gs_name)
    return(gene_sets)
}

#' Run ORA on expression or splicing results
#' @noRd
run_ora <- function(
        paired_diff_result, gene_sets, type, cutoff, min_size){
    
    if(type == "splicing") {
        paired_diff_result <- paired_diff_result[
            !is.na(paired_diff_result$padj_splicing),]
    } else if(type == "expression"){
        paired_diff_result <- paired_diff_result[
            !is.na(paired_diff_result$padj_expression),]
    } else if(type == "paired"){
        paired_diff_result <- paired_diff_result[
            !is.na(paired_diff_result$padj_expression)
            & !is.na(paired_diff_result$padj_expression),]
    }

    universe <- unique(paired_diff_result$gene)
    # Subset significant genes
    sig_genes <- subset_genes(
        paired_diff_result, type = type, cutoff = cutoff)
    
    # ORA on results
    ora <- fgsea::fora(
        gene_sets, genes = sig_genes$gene, 
        universe = universe, minSize = min_size)
    # Add relative risk and enrichment scores
    ora <- compute_enrichment(
        ora, length(sig_genes$gene), length(universe))

    return(ora)
}

#' Run ORA on expression or splicing results
#' @noRd
run_fcs <- function(
        paired_diff_result, gene_sets, type, min_size){
    
    if(type == "splicing") {
        paired_diff_result <- paired_diff_result[
            !is.na(paired_diff_result$padj_splicing),]
    } else if(type == "expression"){
        paired_diff_result <- paired_diff_result[
            !is.na(paired_diff_result$padj_expression),]
    } else if(type == "paired"){
        paired_diff_result <- paired_diff_result[
            !is.na(paired_diff_result$padj_expression)
            & !is.na(paired_diff_result$padj_splicing),]
        paired_diff_result$pvalue_paired <- apply(paired_diff_result[c("pvalue_splicing", "pvalue_expression")], 1, min)
    }
    
    paired_diff_result <- paired_diff_result[!is.na(paired_diff_result$gene),]
    fcs_stats <- -log10(paired_diff_result[, paste0("pvalue_", type)]+0.06)
    fcs_stats <- setNames(fcs_stats, paired_diff_result$gene)
    
    
    # FCS on results
    fcs <- fgsea::fgseaMultilevel(
        pathways = gene_sets,
        stats = fcs_stats,
        nproc = 10,
        scoreType = "pos",
        eps = 10e-320,
        minSize = min_size
    )
    
    colnames(fcs) <- gsub("ES", "enrichment_score", colnames(fcs))
    
    return(fcs)
}

#' Calculate enrichment scores
#' @noRd
compute_enrichment <- function(ora, n_genes, n_universe){
    
    ora$relative_risk <- 
        (ora$overlap / ora$size) /
        (n_genes / n_universe)
    ora$enrichment_score <- log2(ora$relative_risk + 0.06)
    
    return(ora)
}

#' Subset genes to cutoff
#' @noRd
subset_genes <- function(paired_diff_result, type, cutoff){
    
    if(type == "paired"){
        paired_diff_result$padj_splicing[
            is.na(paired_diff_result$padj_splicing)] <- 1
        sig_genes <- paired_diff_result[
            (paired_diff_result[, "padj_expression"] < cutoff
            | (paired_diff_result[, "padj_splicing"] < cutoff)), ]
    } else{ # Expression or splicing only
        padj_col <- ifelse(
            type == "expression", "padj_expression", "padj_splicing")
        
        sig_genes <- paired_diff_result[S4Vectors::complete.cases(
            paired_diff_result[,c(
                paste0("pvalue_", type), "gene", padj_col)]),]
        sig_genes <- sig_genes[
            sig_genes[[padj_col]] < cutoff, ]
    }
    
    return(sig_genes)
}

#' Join ORA results
#' @noRd
join_results <- function(expression, splicing, paired, type = "ora"){
    
    joined <- merge(
        expression,
        splicing,
        by = c("pathway"),
        suffixes = c("_expression", "_splicing"),
        all = TRUE)
    colnames(paired) <- gsub("$", "_paired", colnames(paired))
    colnames(paired)[1] <- "pathway"
    joined <- merge(
        joined,
        paired,
        by = c("pathway"),
        all = TRUE)
    if (type != "ora") return(joined)
    # Fill NAs with zeros
    joined$relative_risk_expression[
        is.na(joined$relative_risk_expression)] <- 0
    joined$enrichment_score_expression[
        is.na(joined$enrichment_score_expression)] <- log2(0.06)
    joined$relative_risk_splicing[
        is.na(joined$relative_risk_splicing)] <- 0
    joined$enrichment_score_splicing[
        is.na(joined$enrichment_score_splicing)] <- log2(0.06)
    joined$relative_risk_splicing[
        is.na(joined$relative_risk_paired)] <- 0
    joined$enrichment_score_splicing[
        is.na(joined$enrichment_score_paired)] <- log2(0.06)
    # Compute shifts in relative risk and enrichment scores
    joined$relative_risk_shift <- 
        joined$relative_risk_splicing - joined$relative_risk_expression
    joined$enrichment_score_shift <- log2(
        abs(joined$relative_risk_shift)) *
        sign(joined$relative_risk_shift)
    joined <- joined[order(joined$enrichment_score_shift), ]
    
    return(joined)
}