#' Inverted versions of in
#' @noRd
#' @examples
#' 1 %!in% 1:10
`%!in%` <- Negate(`%in%`)

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL



#' Wrapper for missing packages
#'
#' @noRd
check_missing_package <- function(package, repo = "CRAN", git_repo = ""){
  
  if (repo == "CRAN"){
    install_function <- "install.packages('"
  } else if (repo == "github") {
    install_function <- paste0("devtools::install_github('", git_repo, "/")
  } else if (repo == "Bioc"){
    install_function <- "BiocManager::install('"
  } else{
    install_function <- "Unknown repository.. "
  }
  
  if(!requireNamespace(package, quietly = TRUE)){
    stop(paste0("Package ", package," is not installed.\n",
                "Please run: ", install_function, package, "')"))
  }
  requireNamespace(package, quietly = TRUE)
}

#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {dir.create(dir_path)}
}

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata"){
  if("data.frame" %in% class(df_colnames)) df_colnames <- colnames(df_colnames)
  if(!is.null(col_name)){
    if(col_name %!in% df_colnames){
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}

#' Check comparison has the right format
#' @noRd
check_comparison <- function(comparison){
  
  if(class(comparison) == "character" & length(comparison) == 1){
    stopifnot("Comparison must have the format 'baseline_v_case' (e.g., '2v1')." = stringr::str_detect(comparison, "v"))
    stopifnot("Comparison must not contain multiple 'v's, the format should be 'baseline_v_case' (e.g., '2v1')." = stringr::str_count(comparison, "v") == 1)
    comparison <- comparison %>% 
      stringr::str_replace(pattern = "_v_", replacement = "v") %>% 
      stringr::str_replace(pattern = " v ", replacement = "v") %>% 
      stringr::str_split(pattern = "v", simplify = TRUE) %>% 
      as.character()
  }
  stopifnot("Comparison must be a character string (e.g., '2v1') or a list (e.g., c('2', '1')). The format being 'baseline_v_case'." = (class(comparison) == "character") & (length(comparison) == 2))
  return(comparison)
}


#' Store result object
#' @noRd
store_result <- function(object, file, analysis = "results", quiet = FALSE){
  check_make_dir("results/")
  if(!startsWith(file, "results")) file <- paste0("results/", file)
  
  if(endsWith(toupper(file), ".RDS")) saveRDS(object, file = file)
  else if(endsWith(tolower(file), ".rdata")) save(object, file = file)
  else if(endsWith(tolower(file), ".csv")) utils::write.csv(object, file =  file)
  else if(endsWith(tolower(file), ".tsv")) utils::write.table(object, file = file)
  else if(endsWith(tolower(file), ".xlsx")) {check_missing_package("writexl"); writexl::write_xlsx(object, path = file)}
  else{stop("Could not store result object ", file)}
  if(!quiet) message("Stored ", analysis, " in ", file)
}



#' Pre-filter
#' @inheritParams paired_de
#' @noRd
pre_filter <- function(tx_count, threshold = 10){
  if(threshold < 1) return(tx_count)
  # Remove low counts
  keep <- rowSums(tx_count) >= threshold
  tx_count <- tx_count[keep,]
  
  return(tx_count)
}

#' Convert character vector into a stats formula
#' @param vector Character vector of variables to add to model
formularise_vector <- function(vector){
  stopifnot("Names must not contain spaces" = stringr::str_detect(vector, " ", negate = TRUE))
  if(is.null(vector)){
    formula <- NULL
  } else if(length(vector) == 1){
    formula <- stats::formula(paste0("~", vector))
  } else{
    formula <- stats::formula(paste0("~", vector[1], paste("+", vector[2:length(vector)], collapse = " ")))
  }
  return(formula)
}

#' Reduce design formula by removing first variable
#' 
#' This internal function assumes the design formula has the format as imposed by \code{\link[pairedGSEA:formularise_vector]{formularise_vector()}}
#'   and that the defining variable is the first occurring variable in the design.
#' @param formula A design formula to reduce. The first variable will be removed.
#' @param formularise (Default: TRUE) A logical determining if the design should be formularised
reduce_formula <- function(formula, formularise = TRUE){
  
  # Convert design formula to character string
  formula_vector <- as.character(formula)[2] %>% 
    stringr::str_split(" \\+ ", simplify = TRUE)
  
  # If design only contains 1 variable, return ~1
  if(length(formula_vector) == 1) reduced_vector <- "1"
  # Else, remove first variable
  else reduced_vector <- formula_vector[2:length(formula_vector)]
  
  if(formularise) return(formularise_vector(reduced_vector))
  
  return(reduced_vector)
  
}


#' Load MSigDB and convert to names list of gene sets
#' 
#' This internal function is wrapper around \code{\link[msigdbr:msigdbr]{msigdbr()}}. Please see their manual for details on its use.
#'   The function extracts the gene set name and a user-defined gene id type (Defualt: "ensembl_gene").
#'   Please make sure the gene IDs match those from your DE analysis.
#'   This function will format the gene sets such that they can be directly used with \code{\link[pairedGSEA:paired_ora]{paired_ora()}}.
#'   
#' @param gene_id_type (Default: "ensemble_gene") The gene ID type to extract. The IDs should match the gene IDs from your DE analysis.
#' @inheritParams msigdbr::msigdbr
#' @keywords internal
prepare_msigdb <- function(gene_id_type = "ensembl_gene",
                           species = "Homo sapiens", 
                           category = "C5"){
  check_missing_package("msigdbr")
  
  gene_sets <- msigdbr::msigdbr(species = species,
                                category = category)
  # Split dataframe based on gene set names
  gene_sets <- gene_sets %>% 
    base::split(x = .[[gene_id_type]], f = .$gs_name)
  return(gene_sets)
}