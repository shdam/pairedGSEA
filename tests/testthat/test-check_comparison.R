test_that("Comparison is correctly reformatted", {
  expect_equal(check_comparison("2v1"), c("2", "1"))
  expect_equal(check_comparison("2 v 1"), c("2", "1"))
  expect_equal(check_comparison("2_v_1"), c("2", "1"))
})
