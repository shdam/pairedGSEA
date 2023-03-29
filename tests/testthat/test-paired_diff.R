test_that("paired_diff works", {
    data("example_se")
    diff_results <- suppressWarnings(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            ))
    expect_s3_class(
        diff_results , "data.frame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
    })

test_that("paired_diff w/limma works", {
    data("example_se")
    diff_results <- suppressWarnings(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            use_limma = TRUE,
            store_results = FALSE,
            quiet = TRUE
        ))
    expect_s3_class(
        diff_results , "data.frame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
})

test_that("Error checks in paired_diff", {
    expect_error(
        paired_diff(
            object = "A string",
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    expect_error(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 4,
            store_results = FALSE,
            quiet = TRUE
            )
        )
    expect_error(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            store_results = TRUE,
            quiet = TRUE
            )
        )
    expect_error(
        paired_diff(
            object = matrix(1:4,2),
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    expect_error(
        paired_diff(
            object = example_se,
            metadata = "metadata/1_GSE154968.xlsx",
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    rownames(example_se) <- 1:nrow(example_se)
    expect_error(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    })

test_that("convert_matrix_to_dds returns DESeqDataSet", {
    counts <- matrix(
        c(3, 0, 4, 5, 0, 7, 0, 2, 4, 1, 0, 5),
        nrow = 4,
        dimnames = list(c("Gene1", "Gene2", "Gene3", "Gene4"), 
                        c("Sample1", "Sample2", "Sample3"))
    )
    metadata <- data.frame(
        row.names = colnames(counts),
        condition = c("Case", "Case", "Baseline")
    )
    design <- stats::model.matrix(~condition, data = metadata)
    dds <- convert_matrix_to_dds(counts, metadata, design)
    expect_s4_class(dds, "DESeqDataSet")
    expect_equal(length(unique(metadata$condition)), 
                 length(unique(colData(dds)$condition)))
})

test_that("prepare_metadata removes irrelevant groups", {
    metadata <- data.frame(
        sample_id = c("A", "B", "C", "D", "E"),
        group = c("irrelevant", "baseline", "case", "case", "case"))
    baseline_case <- c("baseline", "case")
    metadata_filtered <- prepare_metadata(metadata, "group", baseline_case)
    expect_equal(metadata_filtered$group, as.factor(c("baseline", "case", "case", "case")))
})
