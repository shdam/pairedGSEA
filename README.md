
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pairedGSEA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`pariedGSEA` is an R package that helps you to run a paired differential
gene expression (DGE) and splicing (DGS) analysis. Providing a count
matrix, `pariedGSEA` combines the results of `DESeq2` (DGE) and `DEXSeq`
(DGS), aggregates the p-values to gene level, and allows you to run a
subsequent gene set over-representation analysis using its
implementation of the `fgsea::fora` function.

## Installation

``` r
# Install Bioconductor packages with:
# devtools::install_bioc(c("DESeq2", "DEXSeq", "fgsea", "SummarizedExperiment", "S4Vectors", "sva", "BiocParallel"))

# Install pairedGSEA
devtools::install_github("shdam/pairedGSEA")
```

<!-- ## Article -->

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
library("pairedGSEA")
suppressPackageStartupMessages(library("SummarizedExperiment"))

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
#> Number of significant surrogate variables is:  2 
#> Iteration (out of 5 ):1  2  3  4  5
#> converting counts to integer mode
#> Warning in DESeqDataSet(rse, design, ignoreRank = TRUE): some variables in
#> design formula are characters, converting to factors
#> Aggregating p values
diff_results
#> # A tibble: 954 × 7
#>    gene            lfc_deseq pvalue_deseq lfc_dexseq pvalue_de…¹ padj_…² padj_…³
#>    <chr>               <dbl>        <dbl>      <dbl>       <dbl>   <dbl>   <dbl>
#>  1 ENSG00000000419   0.0530   0.483           0.251     3.07e- 1 6.16e-1 4.83e-1
#>  2 ENSG00000001167  -0.383    0.000000317    -0.0478    2.64e- 1 1.88e-6 4.35e-1
#>  3 ENSG00000003402   0.234    0.00000523      0.327     1.28e-10 2.67e-5 2.22e-9
#>  4 ENSG00000005007   0.112    0.00353        -0.689     3.88e- 6 1.02e-2 3.13e-5
#>  5 ENSG00000005020  -0.286    0.00155         0.0996    6.31e- 1 4.94e-3 7.94e-1
#>  6 ENSG00000005810   0.151    0.0235         -0.445     7.66e- 2 5.32e-2 1.61e-1
#>  7 ENSG00000005812   0.319    0.0000132       0.357     2.62e- 1 6.36e-5 4.35e-1
#>  8 ENSG00000005893   0.00183  0.944           0.433     9.67e- 1 9.65e-1 9.90e-1
#>  9 ENSG00000006607   0.102    0.0387          0.302     4.73e-10 8.00e-2 7.00e-9
#> 10 ENSG00000007174   0.463    0.551          -1.05      8.09e- 1 6.76e-1 9.13e-1
#> # … with 944 more rows, and abbreviated variable names ¹​pvalue_dexseq,
#> #   ²​padj_deseq, ³​padj_dexseq
#> # ℹ Use `print(n = ...)` to see more rows
```

Over-representation analysis of results

``` r
# Define gene sets in your preferred way
gene_sets <- pairedGSEA::prepare_msigdb(species = "Homo sapiens", category = "C5", gene_id_type = "ensembl_gene")

ora <- paired_ora(
  paired_diff_result = diff_results,
  gene_sets = gene_sets
  )
#> Identifying differentially expressed genes
#> Running over-representation analysis
#> Joining result
ora
#> # A tibble: 559 × 18
#>    pathway       pval_…¹ padj_…² overl…³ size_…⁴ overl…⁵ size_…⁶ size_…⁷ relat…⁸
#>    <chr>           <dbl>   <dbl>   <int>   <int> <list>    <int>   <int>   <dbl>
#>  1 HP_ONSET        0.670   0.716      20      50 <chr>       401     954   0.952
#>  2 GOCC_GOLGI_M…   0.600   0.664      12      29 <chr>       401     954   0.984
#>  3 GOBP_SMALL_G…   0.538   0.611      12      28 <chr>       401     954   1.02 
#>  4 HP_PTOSIS       0.538   0.611      12      28 <chr>       401     954   1.02 
#>  5 HP_ABNORMALI…   0.472   0.556      12      27 <chr>       401     954   1.06 
#>  6 GOBP_TUBE_DE…   0.374   0.481      20      44 <chr>       401     954   1.08 
#>  7 GOCC_CENTROS…   0.374   0.481      20      44 <chr>       401     954   1.08 
#>  8 HP_HIGH_PALA…   0.332   0.442      16      34 <chr>       401     954   1.12 
#>  9 GOBP_HISTONE…   0.278   0.391      16      33 <chr>       401     954   1.15 
#> 10 HP_ABNORMALI…   0.278   0.391      16      33 <chr>       401     954   1.15 
#> # … with 549 more rows, 9 more variables: enrichment_score_deseq <dbl>,
#> #   pval_dexseq <dbl>, padj_dexseq <dbl>, overlap_dexseq <int>,
#> #   overlapGenes_dexseq <list>, size_genes_dexseq <int>,
#> #   relative_risk_dexseq <dbl>, enrichment_score_dexseq <dbl>,
#> #   enrichment_score_shift <dbl>, and abbreviated variable names ¹​pval_deseq,
#> #   ²​padj_deseq, ³​overlap_deseq, ⁴​size_geneset, ⁵​overlapGenes_deseq,
#> #   ⁶​size_genes_deseq, ⁷​size_universe, ⁸​relative_risk_deseq
#> # ℹ Use `print(n = ...)` to see more rows, and `colnames()` to see all variable names
```

You can now plot the enrichment scores against each other and identify
pathways of interest.

``` r
plot_ora(ora)
```

<img src="man/figures/README-plot-1.png" width="100%" />

## Report issues

If you have any issues or questions regarding the use of pairedGSEA,
please do not hesitate to raise an issue on GitHub. In this way, others
may also benefit from the answers and discussions.
