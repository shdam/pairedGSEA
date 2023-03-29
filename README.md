
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pairedGSEA

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/shdam/pairedGSEA/branch/master/graph/badge.svg)](https://app.codecov.io/gh/shdam/pairedGSEA?branch=master)
<!-- badges: end -->

`pairedGSEA` is an R package that helps you to run a paired differential
gene expression (DGE) and splicing (DGS) analysis. Providing a bulk RNA
count data, `pairedGSEA` combines the results of `DESeq2` (DGE) and
`DEXSeq` (DGS), aggregates the p-values to gene level, and allows you to
run a subsequent gene set over-representation analysis using its
implementation of the `fgsea::fora` function.

## Article

A preprint on `pairedGSEA` is available
[here](https://doi.org/10.1101/2022.08.29.505720).

## Installation

Dependencies

``` r
# Install Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "S4Vectors", "DESeq2", "DEXSeq", "fgsea", "sva", "BiocParallel"))
```

Install `pairedGSEA` from GitHub

``` r
# Install pairedGSEA from github
devtools::install_github("shdam/pairedGSEA", build_vignettes = TRUE)
```

Install `pairedGSEA` from Bioconductor

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("pairedGSEA", version = 'devel')
```

## Documentation

To view documentation for the version of this package installed in your
system, start R and enter:

``` r
browseVignettes("pairedGSEA")
```

## Usage

Please see the User Guide vignette for a detailed description of usage.

Here is a quick runthrough of the functions:

<br> Load example data.

``` r
suppressPackageStartupMessages(library("SummarizedExperiment"))
library("pairedGSEA")

data("example_se")
example_se
#> class: SummarizedExperiment 
#> dim: 5611 6 
#> metadata(0):
#> assays(1): counts
#> rownames(5611): ENSG00000282880:ENST00000635453
#>   ENSG00000282880:ENST00000635195 ... ENSG00000249230:ENST00000504393
#>   ENSG00000249244:ENST00000505994
#> rowData names(0):
#> colnames(6): GSM1499784 GSM1499785 ... GSM1499791 GSM1499792
#> colData names(5): study id source final_description group_nr
```

Run paired differential analysis

``` r
set.seed(500) # For reproducible results

diff_results <- paired_diff(
  example_se,
  group_col = "group_nr",
  sample_col = "id",
  baseline = 1,
  case = 2,
  store_results = FALSE,
  quiet = TRUE
  )
#> Warning in pre_filter(dds, prefilter): 
#> Removing 0 rows with a summed count lower than 10
#> Removing 108 rows with counts in less than 2 samples.
#> No significant surrogate variables
#> converting counts to integer mode
#> Warning in DESeqDataSet(rse, design, ignoreRank = TRUE): some variables in
#> design formula are characters, converting to factors
diff_results
#> # A tibble: 942 × 7
#>    gene            lfc_expression pvalue_expre…¹ padj_…² lfc_s…³ pvalu…⁴ padj_…⁵
#>    <chr>                    <dbl>          <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#>  1 ENSG00000000419        0.0420     0.772       8.46e-1  0.262  2.75e-1 5.12e-1
#>  2 ENSG00000001167       -0.391      0.000000528 3.23e-6 -0.0500 6.83e-1 8.58e-1
#>  3 ENSG00000003402        0.201      0.00000472  2.49e-5  0.338  2.84e-6 3.90e-5
#>  4 ENSG00000005007        0.101      0.0229      5.16e-2  0.484  7.98e-7 1.28e-5
#>  5 ENSG00000005020       -0.301      0.0197      4.49e-2  0.0989 6.02e-1 8.03e-1
#>  6 ENSG00000005810        0.138      0.0405      8.38e-2 -0.524  1.63e-1 3.56e-1
#>  7 ENSG00000005812        0.312      0.0000986   4.19e-4  0.403  3.88e-1 6.41e-1
#>  8 ENSG00000005893       -0.00949    0.751       8.29e-1  0.466  9.51e-1 9.96e-1
#>  9 ENSG00000006607        0.0460     0.0171      3.99e-2  0.334  4.47e-6 5.52e-5
#> 10 ENSG00000007174        0.461      0.339       4.68e-1  0.129  6.91e-1 8.67e-1
#> # … with 932 more rows, and abbreviated variable names ¹​pvalue_expression,
#> #   ²​padj_expression, ³​lfc_splicing, ⁴​pvalue_splicing, ⁵​padj_splicing
```

Over-representation analysis of results

``` r
# Define gene sets in your preferred way
gene_sets <- pairedGSEA::prepare_msigdb(
    species = "Homo sapiens",
    category = "C5",
    gene_id_type = "ensembl_gene"
    )

ora <- paired_ora(
  paired_diff_result = diff_results,
  gene_sets = gene_sets
  )
#> Identifying differentially expressed genes
#> Running over-representation analysis
#> Joining result
ora
#> # A tibble: 555 × 18
#>    pathway       pval_…¹ padj_…² overl…³ size_…⁴ overl…⁵ size_…⁶ size_…⁷ relat…⁸
#>    <chr>           <dbl>   <dbl>   <int>   <int> <list>    <int>   <int>   <dbl>
#>  1 HP_ABNORMAL_…   0.993   0.995      10      38 <chr>       414     942   0.599
#>  2 GOBP_MICROTU…   0.925   0.945       8      25 <chr>       414     942   0.728
#>  3 GOMF_TRANSCR…   0.965   0.974      28      80 <chr>       414     942   0.796
#>  4 HP_ABNORMALI…   0.824   0.858      10      27 <chr>       414     942   0.843
#>  5 HP_ABNORMALI…   0.779   0.825      10      26 <chr>       414     942   0.875
#>  6 HP_ABNORMALI…   0.757   0.813      18      45 <chr>       414     942   0.910
#>  7 GOMF_SEQUENC…   0.722   0.790      34      82 <chr>       414     942   0.943
#>  8 GOBP_REGULAT…   0.637   0.738      14      33 <chr>       414     942   0.965
#>  9 HP_ABNORMAL_…   0.584   0.708      20      46 <chr>       414     942   0.989
#> 10 GOBP_MULTICE…   0.578   0.704      14      32 <chr>       414     942   0.995
#> # … with 545 more rows, 9 more variables: enrichment_score_expression <dbl>,
#> #   pval_splicing <dbl>, padj_splicing <dbl>, overlap_splicing <int>,
#> #   overlapGenes_splicing <list>, size_genes_splicing <int>,
#> #   relative_risk_splicing <dbl>, enrichment_score_splicing <dbl>,
#> #   enrichment_score_shift <dbl>, and abbreviated variable names
#> #   ¹​pval_expression, ²​padj_expression, ³​overlap_expression, ⁴​size_geneset,
#> #   ⁵​overlapGenes_expression, ⁶​size_genes_expression, ⁷​size_universe, …
```

You can now plot the enrichment scores against each other and identify
pathways of interest.

``` r
plot_ora(ora) +
    ggplot2::theme_classic()
```

<img src="man/figures/README-plot-1.png" width="100%" />

## Report issues

If you have any issues or questions regarding the use of pairedGSEA,
please do not hesitate to raise an issue on GitHub. In this way, others
may also benefit from the answers and discussions.
