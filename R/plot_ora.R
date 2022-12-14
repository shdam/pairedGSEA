#' Scatter plot of Over-Representation Analysis results
#' 
#' @param ora Output of \code{\link[pairedGSEA:paired_ora]{paired_ora()}}
#' @param plotly (Default: FALSE) Logical on whether to return plot as an interactive plotly plot or a simple ggplot.
#' @param pattern Highlight pathways containing a specific regex pattern
#' @param cutoff (Default: 0.2) Adjusted p-value cutoff for pathways to include
#' @param lines (Default: TRUE) Whether to show dashed lines
#' @param colors (Default: c("darkgray", "purple", "navy")) Colors to use in plot. The colors are ordered as "Both", "DGS", and "DGE"
#' @family plotting
#' @export
#' @examples 
#' data(example_ora_results)
#' 
#' plot_ora(example_ora_results, pattern = "Telomer")
plot_ora <- function(ora,
                     plotly = FALSE,
                     pattern = NULL,
                     cutoff = 0.05,
                     lines = TRUE,
                     colors = c("darkgray", "purple", "lightblue")){
  
  stopifnot("plot_ora currently only works when Differential Splicing has been run" = any(stringr::str_detect(colnames(ora), "_dexseq")))
  
  check_missing_package(package = "ggplot2")
  if(plotly) check_missing_package(package = "plotly")
  
  # Set default match to FALSE for when pattern is NULL
  ora$pattern_match <- FALSE
  
  # Compute spearman correlation between enrichment scores
  correlation <- ora %>% 
    dplyr::filter(padj_deseq < cutoff | padj_dexseq < cutoff) %>% 
    dplyr::summarise(correlation = stats::cor(enrichment_score_dexseq,
                                              enrichment_score_deseq,
                                              method = "spearman")) %>% 
    dplyr::pull(correlation) %>% 
    round(2)
  
  
  # Extract plot data allowing two layers of dots, one for patter non-matches and one for matches
  plt_data <- ora %>% 
    dplyr::filter(padj_deseq < cutoff | padj_dexseq < cutoff) %>% 
    purrr::when(is.null(pattern) ~ ., # Search for pattern
                TRUE ~ dplyr::mutate(., pattern_match = factor(stringr::str_detect(tolower(pathway), tolower(pattern)), levels = c(FALSE, TRUE)))) %>% 
    dplyr::mutate( # Color dots based on which analysis found the pathway enriched
      plot_color = dplyr::case_when(
        padj_dexseq < 0.05 & padj_deseq < 0.05 ~ "Both",
        padj_dexseq < 0.05 ~ "Only Splicing",
        padj_deseq < 0.05 ~ "Only Expression",
        TRUE ~ "NA"
        ),
      plot_color = factor(plot_color, levels = c("Both", "Only Splicing", "Only Expression"))
    ) %>% 
    dplyr::filter(plot_color != "NA")
  
  # Check that there is data to plot
  stopifnot("No significant gene sets in input data" = nrow(plt_data) > 0)

  # Extract pattern matches (NULL if pattern was not given)
  matches <- plt_data %>% 
    dplyr::filter(pattern_match == TRUE)
  
  if(is.null(pattern)) pattern <- "No pattern"
  
  # Subset colors if not all are needed
  colors <- colors[c("Both", "Only Splicing", "Only Expression") %in% plt_data$plot_color]
  
  # Create plot
  plt <- plt_data %>% 
    # Plot only non-matches
    # dplyr::filter(!pattern_match == FALSE) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = enrichment_score_dexseq,
                 y = enrichment_score_deseq,
                 fill = plot_color,
                 color = plot_color,
                 text = pathway) +
    # ggplot2::geom_point(alpha = 0.3, size = 1.5) +
    # Add correlation text to plot
    ggplot2::annotate("text", label = paste("Spearman's \u03C1:", correlation), x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.3) +
    ggplot2::labs(x = "Gene-Set Enrichment Score\nDifferential Splicing",
                  y = "Gene-Set Enrichment Score\nDifferential Expression",
                  fill = paste("padj <", cutoff),
                  colour = paste("padj <", cutoff)) +
    # Make legend a bit prettier
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3,
                                                                       alpha = 1,
                                                                       shape = c(4, 4, 4)[1:length(colors)],
                                                                       color = colors))) +
    ggplot2::scale_fill_manual(values = colors) +
    # ggplot2::scale_shape_manual(values = 4) +
    ggplot2::scale_colour_manual(values = colors)
  
  if(!(pattern == "No pattern") & nrow(matches) > 0){
    # Add matches to plot if a pattern was given
    plt <- plt +
      ggplot2::geom_point(alpha = 0.3, size = 1.5) +
      ggplot2::aes(text = paste0(pathway, "\nMatch: ", pattern_match),
                   shape = pattern_match) + 
      ggplot2::geom_point(data = matches, size = 2.5, alpha = 0.9, fill = "red", ggplot2::aes(shape = TRUE)) +
      ggplot2::labs(shape = pattern) +
      ggplot2::scale_shape_manual(values = c(4, 23))
    #+ , color = "red"
      # ggplot2::geom_point(data = matches, size = 2.5, alpha = 0.9, ggplot2::aes(shape = FALSE))
  } else{
    plt <- plt +
      ggplot2::geom_point(alpha = 0.3, size = 1.5, shape = 4)
  }
  
  
  if(lines){
    # Plot guide lines 
    plt <- plt +
      ggplot2::geom_abline(slope = 1, color = "gray", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, color = "gray", linetype = "dashed")
  }
  
  
  
  # Make plot interactive if user desires
  if(plotly) plt <- plotly::ggplotly(plt, tooltip = "text")
  
  return(plt)
}