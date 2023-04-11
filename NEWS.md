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
