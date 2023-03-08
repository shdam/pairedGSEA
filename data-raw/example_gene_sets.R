## code to prepare `example_gene_sets` dataset goes here

# Get full list of genesets
gene_sets <- pairedGSEA::prepare_msigdb(
    species = "Homo sapiens",
    category = "C5",
    gene_id_type = "ensembl_gene"
    )

# Subset to telomere genes
set.seed(890)
example_gene_sets <- c(gene_sets[
    stringr::str_detect(tolower(names(gene_sets)),
                        "telomer")],
    sample(gene_sets, 35))


usethis::use_data(example_gene_sets, overwrite = TRUE)
