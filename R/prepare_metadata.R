#' Prepare metadata
#' 
#' @inheritParams paired_gsea
prepare_metadata <- function(metadata, group_col, comparison){
  
  if(typeof(md[1]) == "character" & length(metadata) == 1){
    if(stringr::str_ends(metadata, ".xlsx")) metadata <- readxl::read_excel(metadata)
    else if(stringr::str_ends(metadata, ".csv")) metadata <- readr::read_csv(metadata)
  } else if("data.frame" %!in% class(metadata)){stop("Please provide path to a metadata file or a data.frame object.")}
  
  if(group_col %!in% colnames(metadata)) stop("Could not find column ", group_col, " in metadata.")
  # Ensure comparison is on the right format
  comparison <- pairedGSEA:::check_comparison(comparison)
  # Remove irrelevant groups
  metadata <- metadata[metadata[[group_col]] %in% comparison, ]
  
  # Check metadata content
  if(nrow(metadata) == 0) stop("The values in", group_col, "does not correspond to the comparison values.\nPlease ensure what comes before and after the 'v' in comparison are what is found in the", group_col, "column.")
  
  # Add comparison levels to metadata
  metadata[[group_col]] <- factor(metadata[[group_col]], levels = comparison)
  
  return(metadata)
}