
data("example_se")
test_se <- example_se[1:20, ]
test_dds <- DESeqDataSetFromMatrix(countData = assay(test_se), colData = colData(test_se), design = ~group_nr)

test_that("paired_diff works", {
    diff_results <- suppressWarnings(
        paired_diff(
            object = test_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            experiment_title = "Example",
            store_results = TRUE,
            quiet = FALSE
            ))
    expect_s4_class(
        diff_results , "DFrame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
    expect_true(file.exists("results/Example_dds.RDS"))
    expect_true(file.exists("results/Example_dxd.RDS"))
    expect_true(file.exists("results/Example_expression_results.RDS"))
    expect_true(file.exists("results/Example_splicing_results.RDS"))
    expect_true(file.exists("results/Example_aggregated_pvals.RDS"))
    file.remove("results/Example_dds.RDS")
    file.remove("results/Example_expression_results.RDS")
    file.remove("results/Example_aggregated_pvals.RDS")
    file.remove("results/Example_dxd.RDS")
    file.remove("results/Example_splicing_results.RDS")
    expect_false(file.exists("results/Example_dds.RDS"))
    expect_false(file.exists("results/Example_expression_results.RDS"))
    expect_false(file.exists("results/Example_aggregated_pvals.RDS"))
    expect_false(file.exists("results/Example_dxd.RDS"))
    expect_false(file.exists("results/Example_splicing_results.RDS"))
    })

test_that("paired_diff works with matrix", {
    diff_results <- suppressWarnings(
        paired_diff(
            object = SummarizedExperiment::assay(test_se),
            metadata = SummarizedExperiment::colData(test_se),
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
        ))
    expect_s4_class(
        diff_results , "DFrame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
})

test_that("paired_diff works with rownames", {
    diff_results <- suppressWarnings(
        paired_diff(
            object = test_se,
            group_col = "group_nr",
            sample_col = "rownames",
            baseline = "1",
            case = "2",
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
        ))
    expect_s4_class(
        diff_results , "DFrame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
})

test_that("paired_diff w/limma works", {
    diff_results <- suppressWarnings(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            experiment_title = "Example",
            use_limma = TRUE,
            store_results = FALSE,
            quiet = TRUE
        ))
    expect_s4_class(
        diff_results , "DFrame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
})

test_that("A custom design can be used", {
    test_dds$source <- test_dds$group_nr
    diff_results <- suppressWarnings(
        paired_diff(
            object = test_dds,
            group_col = "source",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            custom_design = stats::formula("~source"),
            quiet = TRUE
        )
    )
    expect_s4_class(
        diff_results , "DFrame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
    
    DESeq2::design(test_dds) <- stats::formula("~source")
    diff_results <- suppressWarnings(
        paired_diff(
            object = test_dds,
            group_col = "source",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            custom_design = TRUE,
            quiet = TRUE
        )
    )
    expect_s4_class(
        diff_results , "DFrame")
    expect_true(
        all(c("padj_splicing", "padj_expression") %in% colnames(diff_results)))
    expect_gt(
        nrow(diff_results), 0)
})

test_that("Error checks in paired_diff", {
    expect_error(
        paired_diff(
            object = test_dds,
            group_col = "group_nr",
            sample_col = "id",
            baseline = 1,
            case = 2,
            experiment_title = "Example",
            store_results = FALSE,
            custom_design = stats::formula("~id"),
            quiet = TRUE
        )
    )
    
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
            object = test_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "22",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    expect_error(
        paired_diff(
            object = test_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
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
            object = test_se,
            metadata = "metadata/1_GSE154968.xlsx",
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    rownames(test_se) <- 1:nrow(test_se)
    expect_error(
        paired_diff(
            object = test_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            experiment_title = "Example",
            store_results = FALSE,
            quiet = TRUE
            )
        )
    })

test_that("User is informed of wrong characters in sample name",{
    test_se$id <- gsub(".$", "-", test_se$id)
    colnames(test_se) <- gsub(".$", "-", colnames(test_se))
    expect_error(
        paired_diff(
            object = test_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
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
    
   # prepare_metadata(metadata, "sample_id", baseline_case)
})

test_that("prepare_metadata errors", {
    expect_error(
        prepare_metadata(
            matrix(1:10), 
            "1",
            "1_2"
    ), "Please provide path to a metadata file or a data.frame /
            DataFrame object.")
    
})

test_data <- data.frame(
    gene = c("A", "A", "B", "B", "C", "C"),
    pvalue = c(0.05, 0.03, 0.01, 0.02, 0.04, 0.06),
    baseMean = c(10, 12, 8, 9, 15, 17),
    lfc = c(1.5, -1.2, 0.8, -0.9, 2.3, -2.1)
)

test_that("aggregate_pvalue handles correct input and produces correct output", {
    # Test default settings
    res <- aggregate_pvalue(test_data)
    expect_s4_class(res, "DFrame")
    expect_equal(nrow(res), 3)
    expect_true(all(c("gene", "lfc", "pvalue") %in% colnames(res)))
    
    # Test custom column names
    test_data_custom <- test_data
    colnames(test_data_custom) <- c("custom_gene", "custom_pvalue", "custom_baseMean", "custom_lfc")
    expect_error(aggregate_pvalue(test_data_custom))
    
    # Test different types
    res_expression <- aggregate_pvalue(test_data, type = "expression")
    res_splicing <- aggregate_pvalue(test_data, type = "splicing")
    expect_false(identical(res_expression, res_splicing))
})

test_that("aggregate_pvalue handles edge cases", {
    # Test empty input
    expect_error(aggregate_pvalue(NULL), "Input data is not a data.frame")
    
    # Test input with missing required columns
    incomplete_data <- test_data[, -1]
    expect_error(aggregate_pvalue(incomplete_data))
    
})



## run_limma

test_that("run_limma handles correct input and produces correct output", {
    res <- run_limma(test_dds, group_col = "group_nr", baseline = "1", case = "2", quiet = TRUE)
    
    expect_type(res, "list")
    expect_equal(length(res), 2)
    expect_true(all(c("expression", "splicing") %in% names(res)))
    
    expect_s4_class(res$expression, "DFrame")
    expect_true(all(c("gene", "transcript", "lfc", "padj", "pvalue", "baseMean") %in% colnames(res$expression)))
    
    expect_s4_class(res$splicing, "DFrame")
    expect_true(all(c("gene", "transcript", "lfc", "padj", "pvalue", "baseMean") %in% colnames(res$splicing)))
})

test_that("run_limma handles result storage", {
    # Test store_results option
    res <- run_limma(test_dds, group_col = "group_nr", baseline = "1", case = "2", store_results = TRUE, quiet = TRUE)
    expect_true(file.exists("results/Experiment Title_limma_results.RDS"))
    file.remove("results/Experiment Title_limma_results.RDS")
})

## run_dexseq

test_that("run_dexseq returns a dataframe with correct column names", {
    result <- suppressWarnings(
        run_dexseq(
            test_dds,
            group_col = "group_nr",
            baseline = "1",
            case = "2",
            experiment_title = "test",
            store_results = FALSE,
            quiet = TRUE,
            parallel = FALSE))
    
    expect_s4_class(result, "DFrame")
    expected_colnames <-  c(
        "gene", "transcript", "baseMean", "lfc", "padj", "pvalue")
    expect_equal(colnames(result), expected_colnames)
})

test_that("run_dexseq raises an error if rownames are not in gene:transcript format", {
    # Modify rownames to break the format
    wrong_rownames <- gsub(":", "_", rownames(test_dds))
    rownames(test_dds) <- wrong_rownames
    
    expect_error(run_dexseq(
        test_dds,
        group_col = "group_nr",
        baseline = "1",
        case = "2",
        experiment_title = "test",
        store_results = FALSE,
        quiet = TRUE,
        parallel = FALSE),
        "Please ensure the rownames have the format 'gene:transcript'")
})


## run_deseq
expression_results <- run_deseq(
    test_dds,
    group_col = "group_nr",
    baseline = "1",
    case = "2",
    quiet = TRUE
)
test_that("run_deseq returns the correct structure", {
    
    expect_s4_class(expression_results, "DFrame")
    expect_named(expression_results, c(
        "gene", "transcript", "lfc", "pvalue", "padj", "baseMean"))
})

test_that("run_deseq returns correct values", {

    expect_type(expression_results$gene, "character")
    expect_type(expression_results$transcript, "character")
    expect_type(expression_results$lfc, "double")
    expect_type(expression_results$pvalue, "double")
    expect_type(expression_results$padj, "double")
    expect_type(expression_results$baseMean, "double")
    
    expect_true(all(!is.na(expression_results$gene)))
    expect_true(all(!is.na(expression_results$transcript)))
})

test_that("run_deseq handles invalid arguments", {
    expect_error(run_deseq(
        test_dds, group_col = "invalid", baseline = "1", case = "2", quiet = TRUE),
        "invalid should be the name of a factor in the colData of the DESeqDataSet")
})
