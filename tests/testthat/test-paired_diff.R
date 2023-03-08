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
    all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
  expect_gt(
    nrow(diff_results), 0)
})

test_that("Error checks in paired_diff",{
  expect_error(
    paired_diff("A string",
                group_col = "group_nr",
                sample_col = "id",
                baseline = 1,
                case = 2,
                experiment_title = "Example",
                store_results = FALSE)
  )
    expect_error(
        paired_diff(example_se,
                    group_col = "group_nr",
                    sample_col = "id",
                    baseline = 1,
                    case = 4,
                    store_results = FALSE)
    )
  expect_error(
    paired_diff(example_se,
                group_col = "group_nr",
                sample_col = "id",
                baseline = 1,
                case = 2,
                store_results = TRUE)
  )
  expect_error(
    paired_diff(matrix(1:4,2),
                group_col = "group_nr",
                sample_col = "id",
                baseline = 1,
                case = 2,
                experiment_title = "Example",
                store_results = FALSE)
  )
  expect_error(
    paired_diff(example_se,
                metadata = "metadata/1_GSE154968.xlsx",
                group_col = "group_nr",
                sample_col = "id",
                baseline = 1,
                case = 2,
                experiment_title = "Example",
                store_results = FALSE)
  )
  rownames(example_se) <- 1:nrow(example_se)
  expect_error(
    paired_diff(example_se,
                group_col = "group_nr",
                sample_col = "id",
                baseline = 1,
                case = 2,
                experiment_title = "Example",
                store_results = FALSE)
  )
})
