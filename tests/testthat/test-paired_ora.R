test_that("paired_ora works", {
    data("example_diff_result")
    data("example_gene_sets")
  
    ora_test <- paired_ora(example_diff_result, example_gene_sets)
    expect_s3_class(
        ora_test, "data.frame"
        )
  
})
