## code to prepare `example_ora_results` dataset goes here


data("example_diff_result")
data("example_gene_sets", package = "pairedGSEA")

example_ora_results <- paired_ora(
    example_diff_result,
    example_gene_sets)

usethis::use_data(example_ora_results, overwrite = TRUE)
