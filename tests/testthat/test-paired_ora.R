test_that("paired_ora works", {
  data("example_diff_result")
  gene_sets <- pairedGSEA::prepare_msigdb(species = "Homo sapiens", category = "C5", gene_id_type = "ensembl_gene")
  
  ora_test <- paired_ora(example_diff_result, example_gene_sets, store_results = FALSE)
  expect_s3_class(
    ora_test, "data.frame"
  )
  
})
