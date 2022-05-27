#' Scatter plot of Over-Representation Analysis results
#' 
#' @param ora Output of \code{\link[pairedGSEA:paired_ora]{paired_ora()}}
#' @param plotly (Default: FALSE) Logical on whether to return plot as an interactive plotly plot or a simple ggplot.
#' @param pattern Highlight pathways containing a specific regex pattern
#' @param cutoff (Default: 0.2) Adjusted p-value cutoff for pathways to include
#' @param lines (Default: TRUE) Whether to show dashed lines
#' @export
plot_ora <- function(ora, plotly = FALSE, pattern = NULL, cutoff = 0.05, lines = TRUE){
  
  check_missing_package(package = "ggplot2")
  if(plotly) check_missing_package(package = "plotly")
  
  # Set default match to FALSE for when pattern is NULL
  ora$pattern_match <- FALSE
  
  plt <- ora %>% 
    dplyr::filter(padj_deseq < cutoff | padj_dexseq < cutoff) %>% 
    purrr::when(is.null(pattern) ~ ., # Search for pattern
                TRUE ~ dplyr::mutate(., pattern_match = factor(stringr::str_detect(tolower(pathway), tolower(pattern)), levels = c(FALSE, TRUE)))) %>% 
    ggplot2::ggplot() +
    ggplot2::aes(x = enrichment_score_dexseq,
                 y = enrichment_score_deseq,
                 color = pattern_match,
                 shape = pattern_match,
                 alpha = pattern_match,
                 text = pathway) +
    ggplot2::geom_point() + 
    ggplot2::coord_cartesian(xlim = c(-4, 6), ylim = c(-4,6)) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::scale_shape_manual(values = c(23, 18)) +
    ggplot2::scale_alpha_manual(values = c(0.1, 0.9)) +
    ggplot2::labs(x = "Enrichment Score - Differential Splicing",
                  y = "Enrichment Score - Differential Expression",
                  color = pattern,
                  shape = pattern,
                  alpha = pattern) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::theme_classic()
  
  if(lines){
    plt <- plt +
      ggplot2::geom_abline(slope = 1, color = "gray", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, color = "gray", linetype = "dashed")
  }
  
  # Remove legend if no pattern
  if(is.null(pattern)) plt <- plt + ggplot2::theme(legend.position = "none")
  if(plotly) plt <- plotly::ggplotly(plt, tooltip = "text")
  
  return(plt)
}