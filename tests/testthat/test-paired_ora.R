test_that("paired_ora works", {
    data("example_diff_result")
    data("example_gene_sets")
  
    ora_test <- paired_ora(example_diff_result, example_gene_sets)
    expect_s3_class(
        ora_test, "data.frame"
        )
  
})


test_that("prepare_msigdb returns a list of gene sets", {
    gene_sets <- prepare_msigdb()
    expect_true(is.list(gene_sets))
    expect_true(length(gene_sets) > 0)
    expect_true(all(sapply(gene_sets, is.character)))
})
