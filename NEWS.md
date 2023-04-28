# pairedGSEA 0.99.6

* Improved test coverage and depth
* Coerced `paired_ora` output to `DataFrame` from `data.table`
* Added value field to data man pages


# pairedGSEA 0.99.5

* Remove non-exported man pages

# pairedGSEA 0.99.4

* Remove filter_gene_sets option, as it is inherent in `fgsea::ora`
* Added test depth and coverage for `paired_diff` and `paired_ora`
* Fixed wrong storing location for splicing intermediate results

# pairedGSEA 0.99.3

* Rewrote code base to remove tidyverse dependencies
* Added `\code{...}` and `\link{...}` where appropriate in documentation
* Added input parameter checks
* Reduced redundant input parameters from `aggregate_pvalue`
* Modularized `paired_ora` and `plot_ora`
* Added `filter_gene_sets` parameter to help users reduce gene set bias
* Increased test coverage and depth
* Moved data scripts to `inst/script`

# pairedGSEA 0.99.2

* Implement `limma` as alternative analysis method
* Increased test coverage significantly
* Removed usage of deprecated `purrr::when`

# pairedGSEA 0.99.1

* Updated vignette with a brief paragraph on the motivation of `pairedGSEA`.
* Updated NAMESPACE to include all imported packages.
Suggested packages are added in notes.

# pairedGSEA 0.99.0

* Submitted to Bioconductor
