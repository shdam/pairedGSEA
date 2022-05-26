#' Scatter plot of Over-Representation Analysis results
#' 
#' @param ora Output of \code{\link[pairedGSEA:paired_ora]{paired_ora()}}
#' @param plotly (Default: FALSE) Logical on whether to return plot as an interactive plotly plot or a simple ggplot.
#' @param pattern Highlight pathways containing a specific regex pattern
#' @export
plot_ora <- function(ora, plotly = FALSE, pattern = NULL){
  
  check_missing_package(package = "ggplot2")
  if(plotly) check_missing_package(package = "plotly")
  
  # Set default match to FALSE for when pattern is NULL
  ora$pattern_match <- FALSE
  
  plt <- ora %>% 
    dplyr::filter(padj_deseq < 0.05 | padj_dexseq < 0.05) %>% 
    purrr::when(is.null(pattern) ~ ., # Search for pattern
                TRUE ~ dplyr::mutate(., pattern_match = factor(stringr::str_detect(tolower(pathway), tolower(pattern)), levels = c(FALSE, TRUE)))) %>% 
    ggplot2::ggplot() +
    ggplot2::aes(x = log2(relative_risk_dexseq+0.06),
                 y = log2(relative_risk_deseq+0.06),
                 color = pattern_match,
                 text = pathway) +
    ggplot2::geom_point() + 
    ggplot2::geom_abline(slope = 1, color = "gray", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = c(-4, 6), ylim = c(-4,6)) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::labs(x = "Enrichment Score - Differential Splicing",
                  y = "Enrichment Score - Differential Expression",
                  color = pattern) +
    ggplot2::theme_classic()
  
  # Remove legend if no pattern
  if(is.null(pattern)) plt <- plt + ggplot2::theme(legend.position = "none")
  if(plotly) plt <- plotly::ggplotly(plt, tooltip = "text")
  
  return(plt)
}