data("example_ora_results")
test_that("plot_ora generates a ggplot object", {
    plt <- plot_ora(example_ora_results)
    expect_s3_class(plt, "ggplot")
})

test_that("plot_ora works with a DFrame", {
    plt <- plot_ora(S4Vectors::DataFrame(example_ora_results))
    expect_s3_class(plt, "ggplot")
})

test_that("plot_ora returns an interactive plotly plot when plotly = TRUE", {
    plt <- plot_ora(example_ora_results, plotly = TRUE)
    expect_s3_class(plt, "plotly")
})

test_that("plot_ora throws an error if no significant gene sets in input data", {
    ora_empty <- example_ora_results[0, ]
    expect_error(plot_ora(ora_empty),
                 "No over-represented gene sets found.")
})

test_that("plot_ora throws warning if pattern not found", {
    expect_warning(plot_ora(example_ora_results, pattern = 123),
                 "No matches found with pattern: 123")
})

test_that("plot_ora throws an error if input data is missing splicing data",
          {
              ora_missing_splicing <- example_ora_results[, -grep(
                  "_splicing", colnames(example_ora_results))]
              expect_error(
                  plot_ora(ora_missing_splicing),
                  "plot_ora currently only works when Differential Splicing has been run"
              )
          })

test_that("plot_ora highlights pattern matches correctly",
          {
              plt <- plot_ora(example_ora_results, pattern = "TELOMERE")
              # Check that at least one match was found
              expected_matches <- grepl("TELOMERE", plt$data$pathway)
              expect_true(any(expected_matches))
              # Check that only matches are red triangles
              matches <- plt$data$pattern_match
              # Count the number of matches and non-matches
              match_counts <- table(matches)
              # Check that the expected numbers of matches and non-matches are found
              expect_equal(as.list(match_counts),
                           list(
                               "FALSE" = sum(!expected_matches),
                               "TRUE" = sum(expected_matches)))
          })


test_that("plot_ora returns a plot with guide lines when lines = TRUE",
          {
              plt <- plot_ora(example_ora_results, lines = TRUE)
              # Check that guide lines are present
              expect_equal(length(plt$layers), 5)
              expect_equal("slope", plt$layers[[3]]$geom$required_aes[1])
              expect_equal("yintercept", plt$layers[[4]]$geom$required_aes)
              expect_equal("xintercept", plt$layers[[5]]$geom$required_aes)
          })

test_that("plot_ora returns a plot without guide lines when lines = FALSE",
          {
              plt <- plot_ora(example_ora_results, lines = FALSE)
              # Check that guide lines are not present
              expect_equal(length(plt$layers), 2)
          })

