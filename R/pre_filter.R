#' Pre-filter
#' @inheritParams paired_gsea
#' @noRd
pre_filter <- function(tx_count, threshold = 10){
  if(threshold < 1) return(tx_count)
  # Remove low counts
  keep <- rowSums(tx_count) >= threshold
  tx_count <- tx_count[keep,]
  
  return(tx_count)
}