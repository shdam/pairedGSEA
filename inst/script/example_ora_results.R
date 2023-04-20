## code to prepare `example_ora_results` dataset goes here


data("example_diff_result")
gene_sets <- pairedGSEA::prepare_msigdb(
    species = "Homo sapiens", 
    category = "C5", 
    gene_id_type = "ensembl_gene"
)

example_ora_results <- paired_ora(
    example_diff_result,
    gene_sets)

usethis::use_data(example_ora_results, overwrite = TRUE)


