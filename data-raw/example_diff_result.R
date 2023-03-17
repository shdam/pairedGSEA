## code to prepare `example_diff_result` dataset goes here

library("pairedGSEA")
set.seed(500) # For reproducible results
data("example_se", package = "pairedGSEA")
example_diff_result <- diff_results <- paired_diff(
  object = example_se,
  group_col = "group_nr",
  sample_col = "id",
  baseline = 1,
  case = 2,
  store_results = FALSE
)

usethis::use_data(example_diff_result, overwrite = TRUE)
