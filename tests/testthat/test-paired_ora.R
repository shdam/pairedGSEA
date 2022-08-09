test_that("paired_ora works", {
  data("example_diff_result")
  gene_sets <- pairedGSEA::prepare_msigdb()
  
  ora_test <- paired_ora(example_diff_result, gene_sets)
  expect_s3_class(
    ora_test, "data.frame"
  )
  
})
