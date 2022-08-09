test_that("paired_diff works", {
  data(example_se)
  diff_results <- suppressWarnings(
    paired_diff(
    object = example_se,
    group_col = "group_nr",
    sample_col = "id",
    baseline = 1,
    case = 2,
    experiment_title = "Example",
    store_results = FALSE
  ))
  expect_s3_class(
    diff_results , "data.frame")
  expect_true(
    all(c("padj_dexseq", "padj_deseq") %in% colnames(diff_results)))
  expect_gt(
    nrow(diff_results), 0)
})
