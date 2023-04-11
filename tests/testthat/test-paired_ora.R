data("example_diff_result")
data("example_gene_sets")



test_that("paired_ora works", {
    ora_test <- paired_ora(example_diff_result, example_gene_sets)
    expect_equal(nrow(ora_test), 8)
})

test_that("paired_ora works with limma", {
    data("example_se")
    diff_results <- suppressWarnings(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            use_limma = TRUE
        ))
    
    ora_test <- paired_ora(diff_results, example_gene_sets)
    expect_equal(nrow(ora_test), 8)
})

test_that("prepare_msigdb returns a list of gene sets", {
    gene_sets <- prepare_msigdb()
    expect_true(is.list(gene_sets))
    expect_true(length(gene_sets) > 0)
    expect_true(all(sapply(gene_sets, is.character)))
})

genes <- c("A", "B", "C")
gene_sets <- list(c("A", "D", "E"), c("F", "G", "B"), c("H", "I", "J"))

# Define test cases
test_that("rm_gene_sets removes gene sets without any genes", {
    expect_equal(length(rm_gene_sets(gene_sets, genes)), 2)
})

test_that("rm_gene_sets returns the correct gene sets", {
    expect_equal(rm_gene_sets(gene_sets, genes), list(c("A", "D", "E"), c("F", "G", "B")))
})

test_that("it subsets genes to a cutoff based on type", {
    expression_genes <- subset_genes(example_diff_result, "expression", 0.05)
    expect_true(all(expression_genes$padj_expression < 0.05))
    
    splicing_genes <- subset_genes(example_diff_result, "splicing", 0.05)
    expect_true(all(splicing_genes$padj_splicing < 0.05))
})

# Create a mock ORA data frame
ora <- data.frame(
    pathway = c("Pathway A", "Pathway B", "Pathway C"),
    overlap = c(10, 10, 10),
    size = c(100, 200, 300),
    stringsAsFactors = FALSE
)

# Define test cases
test_that("compute_enrichment returns correct enrichment scores", {
    # Test with n_genes = 1000 and n_universe = 10000
    result <- compute_enrichment(ora, n_genes = 1000, n_universe = 10000)
    expect_equal(result$relative_risk, c(1, 0.5, 1/3))
    expect_equal(result$enrichment_score, c(log2(1+0.06), log2(0.5+0.06), log2(1/3+0.06)))
})
