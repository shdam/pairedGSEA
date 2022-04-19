#' Inverted versions of in
#'
#' @noRd
#' @examples
#' 1 %!in% 1:10
`%!in%` <- Negate(`%in%`)

#' Wrapper for missing packages
#'
#' @noRd
check_missing_package <- function(package, repo = "CRAN", git_repo = ""){
  
  if (repo == "CRAN"){
    install_function <- "install.packages('"
  } else if (repo == "github") {
    install_function <- paste0("devtools::install_github('", git_rep, "/")
  } else if (repo == "Bioc"){
    install_function <- "BiocManager::install('"
  }
  
  if(package %!in% rownames(installed.packages())){
    stop(paste0("Package ", package," is not installed.\n",
                "Please run: ", install_function, package, "')"))
  }
  requireNamespace(package, quietly = TRUE)
}

#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir.path) {
  if (!dir.exists(dir.path)) {dir.create(dir.path)}
}

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata"){
  if(!is.null(col_name)){
    if(col_name %!in% df_colnames){
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}

#' Check comparison has the right format
check_comparison <- function(comparison){
  
  if(class(comparison) == "character" & length(comparison) == 1){
    stopifnot("Comparison must have the format 'baseline_v_case' (e.g., '2v1')." = stringr::str_detect(comparison, "v"))
    stopifnot("Comparison must not contain multiple 'v's, the format should be 'baseline_v_case' (e.g., '2v1')." = stringr::str_count(comparison, "v") == 1)
    comparison <- comparison %>% 
      stringr::str_replace(pattern = "_v_", replacement = "v") %>% 
      stringr::str_split(pattern = "v", simplify = TRUE) %>% 
      as.character()
  }
  stopifnot("Comparison must be a character string (e.g., '2v1') or a list (e.g., c('2', '1')). The format being 'baseline_v_case'." = (class(comparison) == "character") & (length(comparison) == 2))
  return(comparison)
}


#' Store result object
#' @noRd
store_result <- function(object, file, analysis = "results"){
  pairedGSEA:::check_make_dir("results/")
  if(!startsWith(file, "results")) file <- paste0("results/", file)
  
  if(endsWith(toupper(file), ".RDS")) saveRDS(object, file = file)
  else if(endsWith(tolower(file), ".rdata")) save(object, file = file)
  else if(endsWith(tolower(file), ".csv")) write.csv(object, file =  file)
  else if(endsWith(tolower(file), ".tsv")) write.table(object, file = file)
  else if(endsWith(tolower(file), ".xlsx")) {pairedGSEA:::check_missing_package("openxlsx"); openxlsx::write.xlsx(object, file = file)}
  else{stop("Could not store result object ", file)}
  message("Stored ", analysis, " in ", file)
}

#' Convert matrix to DESeqDataSet
#' @inheritParams paired_gsea
convert_matrix_to_dds <- function(tx_count, metadata, group_col){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = tx_count,
                                 colData = metadata,
                                 design = as.formula(paste0("~", group_col)))
  return(dds)
}