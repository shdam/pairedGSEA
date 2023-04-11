# Tests for check_missing_package ----
test_that("check_missing_package throws an error when package is missing",
          {
              expect_error(check_missing_package("nonexistingpackage"))
          })
test_that("check_missing_package does not throw an error when package is installed",
          {
              expect_no_error(check_missing_package("dplyr"))
          })

# Tests for check_make_dir ----
test_that("check_make_dir creates a directory when it doesn't exist",
          {
              dir_path <- "test_dir"
              if (dir.exists(dir_path)) unlink(dir_path, recursive = TRUE)
              expect_false(dir.exists(dir_path))
              check_make_dir(dir_path)
              expect_true(dir.exists(dir_path))
              unlink(dir_path, recursive = TRUE)
          })
test_that("check_make_dir does nothing when directory already exists",
          {
              dir_path <- "test_dir"
              if (dir.exists(dir_path)) unlink(dir_path, recursive = TRUE)
              expect_false(dir.exists(dir_path))
              dir.create(dir_path)
              expect_true(dir.exists(dir_path))
              check_make_dir(dir_path)
              expect_true(dir.exists(dir_path))
              unlink(dir_path, recursive = TRUE)
          })
# Tests for check_colname ----
test_that("check_colname throws an error when column is not in dataframe",
          {
              df_colnames <- c("col1", "col2", "col3")
              col_name <- "nonexistingcolumn"
              expect_error(check_colname(df_colnames, col_name))
          })
test_that("check_colname does not throw an error when column is in dataframe",
          {
              df_colnames <- c("col1", "col2", "col3")
              col_name <- "col2"
              expect_no_error(check_colname(df_colnames, col_name))
          })

# Tests for check_comparison ----
test_that("check_comparison works with valid comparison as character string", {
    comparison <- "1v2"
    comparison2 <- "1 vs 2"
    expect_equal(check_comparison(comparison), c("1", "2"))
    expect_equal(check_comparison(comparison2), c("1", "2"))
})

test_that("check_comparison works with valid comparison as list", {
    comparison <- list("1", "2")
    expect_equal(check_comparison(comparison), as.character(comparison))
})

test_that("check_comparison returns error with invalid comparison as character string", {
    comparison <- "1 s 2"
    expect_error(check_comparison(comparison), 
                 "Comparison must have the format 'baseline_v_case'")
})

test_that("check_comparison returns error with invalid comparison as list", {
    comparison <- list("1", "2", "3")
    expect_error(check_comparison(comparison))
})

test_that("check_comparison returns error with missing comparison", {
    comparison <- NULL
    expect_error(check_comparison(comparison))
})


# Test formularise_vector ----
test_that("formularise_vector converts character vector to formula", {
    vec <- c("treatment", "time", "subject")
    formula <- formularise_vector(vec)
    expect_equal(as.character(formula)[2], "treatment + time + subject")
})

test_that("formularise_vector returns NULL for empty vector", {
    vec <- NULL
    formula <- formularise_vector(vec)
    expect_null(formula)
})

# Test reduce_formula ----
test_that("reduce_formula removes first variable from formula", {
    formula <- as.formula("~ treatment + time + subject")
    reduced_formula <- reduce_formula(formula)
    expect_equal(as.character(reduced_formula)[2], "time + subject")
})

test_that("reduce_formula returns NULL for formula with only one variable", {
    formula <- as.formula("~treatment")
    reduced_formula <- reduce_formula(formula, formularise = FALSE)
    expect_equal(reduced_formula, "1")
})

# Test convert_matrix_to_dds and pre_filter
test_that("convert_matrix_to_dds converts matrix to DESeq2 object", {
    # Create example count matrix
    tx_count <- matrix(c(10, 0, 20, 15, 50, 30), nrow = 3, ncol = 2,
                       dimnames = list(c("ENSG001", "ENSG002", "ENSG003"),
                                       c("sample1", "sample2")))
    
    # Create example metadata
    metadata <- data.frame(row.names = colnames(tx_count),
                           condition = c("control", "case"))
    
    # Create example design formula
    design_formula <- "~condition"
    design <- stats::formula(design_formula)
    
    # Convert matrix to DESeq2 object
    dds <- suppressWarnings(convert_matrix_to_dds(tx_count, metadata, design))
    
    # Check if object is a DESeq2DataSet
    expect_s4_class(dds, "DESeqDataSet")
    # Check if object contains count data
    expect_equal(DESeq2::counts(dds), tx_count)
    # Check if object contains sample metadata
    expect_s4_class(SummarizedExperiment::colData(dds), "DataFrame")
    # Check if object contains design formula
    expect_equal(dds@design, design)
})

# Test pre_filter
test_that("convert_matrix_to_dds converts matrix to DESeq2 object", {
    # Create example count matrix
    tx_count <- matrix(c(10, 0, 20, 15, 50, 30), nrow = 3, ncol = 2,
                       dimnames = list(c("ENSG001", "ENSG002", "ENSG003"),
                                       c("sample1", "sample2")))
    
    # Create example metadata
    metadata <- data.frame(row.names = colnames(tx_count),
                           condition = c("control", "case"))
    
    # Create example design formula
    design_formula <- "~condition"
    design <- as.formula(design_formula)
    
    # Convert matrix to DESeq2 object
    dds <- suppressWarnings(convert_matrix_to_dds(tx_count, metadata, design))
    
    dds_filtered <- suppressWarnings(pre_filter(dds, threshold = 30))
    expect_equal(dim(dds_filtered), c(1, 2))
    dds_filtered <- suppressWarnings(pre_filter(dds, threshold = 1))
    expect_equal(dim(dds_filtered), c(2, 2))
    dds_filtered <- suppressWarnings(pre_filter(dds, threshold = 0))
    expect_equal(dim(dds_filtered), c(3, 2))
})


test_that("store_result function works correctly", {
    # create a dummy object to save
    x <- data.frame(matrix(1:6, ncol = 2))
    colnames(x) <- c("A", "B")
    # test RDS
    store_result(x, "test_result.RDS")
    expect_true(file.exists("results/test_result.RDS"))
    file.remove("results/test_result.RDS")
    expect_false(file.exists("results/test_result.RDS"))
    # test csv
    store_result(x, "results/test_result.csv")
    expect_true(file.exists("results/test_result.csv"))
    file.remove("results/test_result.csv")
    expect_false(file.exists("results/test_result.csv"))
    # test xlsx
    store_result(x, "results/test_result.xlsx")
    expect_true(file.exists("results/test_result.xlsx"))
    file.remove("results/test_result.xlsx")
    expect_false(file.exists("results/test_result.xlsx"))
    # test rdata
    store_result(x, "results/test_result.rdata")
    expect_true(file.exists("results/test_result.rdata"))
    file.remove("results/test_result.rdata")
    expect_false(file.exists("results/test_result.rdata"))
    # test tsv
    store_result(x, "results/test_result.tsv")
    expect_true(file.exists("results/test_result.tsv"))
    file.remove("results/test_result.tsv")
    expect_false(file.exists("results/test_result.tsv"))
})
