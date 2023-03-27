#' Scatter plot of Over-Representation Analysis results
#' 
#' @param ora Output of \code{\link[pairedGSEA:paired_ora]{paired_ora()}}
#' @param plotly (Default: FALSE) Logical on whether to return plot as an
#' interactive plotly plot or a simple ggplot.
#' @param pattern Highlight pathways containing a specific regex pattern
#' @param cutoff (Default: 0.2) Adjusted p-value cutoff for pathways to include
#' @param lines (Default: TRUE) Whether to show dashed lines
#' @param colors (Default: c("darkgray", "purple", "navy"))
#' Colors to use in plot. The colors are ordered as "Both", "DGS", and "DGE"
#' @importFrom stringr str_detect
#' @importFrom dplyr filter summarise pull mutate case_when
#' @importFrom stats cor
#' @importFrom ggplot2 ggplot aes annotate labs guides guide_legend
#' scale_fill_manual scale_colour_manual geom_point scale_shape_manual
#' geom_abline geom_hline geom_vline
#' @note Suggested: importFrom plotly ggplotly
#' @family plotting
#' @export
#' @return A ggplot
#' @usage
#' plot_ora(
#'     ora,
#'     plotly = FALSE,
#'     pattern = NULL,
#'     cutoff = 0.05,
#'     lines = TRUE,
#'     colors = c("darkgray", "purple", "lightblue")
#'     )
#' @examples 
#' data(example_ora_results)
#' 
#' plot_ora(example_ora_results, pattern = "Telomer")
plot_ora <- function(
        ora,
        plotly = FALSE,
        pattern = NULL,
        cutoff = 0.05,
        lines = TRUE,
        colors = c("darkgray", "purple", "lightblue")){

    if(is(pattern, "numeric")) pattern <- as.character(pattern)
    stopifnot(
        "Pattern must be a string"
        = is.null(pattern) || is(pattern, "character"))
    
    stopifnot(
        "plot_ora currently only works when Differential Splicing has been run"
        = any(stringr::str_detect(colnames(ora), "_splicing")))
    
    check_missing_package(package = "ggplot2")
    if(plotly) check_missing_package(package = "plotly")
    
    # Set default match to FALSE for when pattern is NULL
    ora$pattern_match <- FALSE
    
    # Filter ora
    ora <- ora %>% 
    dplyr::filter(padj_expression < cutoff | padj_splicing < cutoff)
    
    # Compute spearman correlation between enrichment scores
    correlation <- ora %>% 
    dplyr::summarise(correlation = stats::cor(
        enrichment_score_splicing,
        enrichment_score_expression,
        method = "spearman")) %>% 
    dplyr::pull(correlation) %>% 
    round(2)
    
    
    # Extract plot data allowing two layers of dots, one for patter non-matches
    # and one for matches
    if(!is.null(pattern)) {
        ora$pattern_match <- factor(
            stringr::str_detect(tolower(ora$pathway), tolower(pattern)),
            levels = c(FALSE, TRUE))
    }
    
    plt_data <- ora %>% 
        dplyr::mutate( 
        # Color dots based on which analysis found the pathway enriched
            plot_color = dplyr::case_when(
                padj_splicing < cutoff & padj_expression < cutoff ~ "Both",
                padj_splicing < cutoff ~ "Only Splicing",
                padj_expression < cutoff ~ "Only Expression",
                TRUE ~ "NA"
                ),
            plot_color = factor(
                plot_color, levels = c(
                    "Both", "Only Splicing",
                    "Only Expression"))
            ) %>% 
        dplyr::filter(plot_color != "NA")
    # Check that there is data to plot
    stopifnot("No significant gene sets in input data" = nrow(plt_data) > 0)
    
    # Extract pattern matches (NULL if pattern was not given)
    matches <- plt_data %>% 
        dplyr::filter(pattern_match == TRUE)
    
    if(is.null(pattern)) pattern <- "No pattern"
    
    # Subset colors if not all are needed
    colors <- colors[
        c("Both", "Only Splicing", "Only Expression") %in%
            plt_data$plot_color]
    
    # Create plot
    plt <- plt_data %>% 
    # Plot only non-matches
    # dplyr::filter(!pattern_match == FALSE) %>%
        ggplot2::ggplot() +
        ggplot2::aes(
            x = enrichment_score_splicing,
            y = enrichment_score_expression,
            fill = plot_color,
            color = plot_color,
            text = pathway) +
    # ggplot2::geom_point(alpha = 0.3, size = 1.5) +
    # Add correlation text to plot
        ggplot2::annotate(
            "text", label = paste("Spearman's \u03C1:", correlation),
            x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.3) +
        ggplot2::labs(
            x = "Gene-Set Enrichment Score\nDifferential Splicing",
            y = "Gene-Set Enrichment Score\nDifferential Expression",
            fill = paste("padj <", cutoff),
            colour = paste("padj <", cutoff)) +
    # Make legend a bit prettier
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(
            size = 3,
            alpha = 1,
            shape = c(4, 4, 4)[seq_len(length(colors))],
            color = colors))) +
        ggplot2::scale_fill_manual(values = colors) +
        # ggplot2::scale_shape_manual(values = 4) +
        ggplot2::scale_colour_manual(values = colors)
    
    if(!(pattern == "No pattern") & nrow(matches) > 0){
    # Add matches to plot if a pattern was given
        plt <- plt +
            ggplot2::geom_point(alpha = 0.3, size = 1.5) +
            ggplot2::aes(
                text = paste0(pathway, "\nMatch: ", pattern_match),
                shape = pattern_match) + 
            ggplot2::geom_point(
                data = matches, size = 2.5, alpha = 0.9, fill = "red",
                ggplot2::aes(shape = TRUE)) +
            ggplot2::labs(shape = pattern) +
            ggplot2::scale_shape_manual(values = c(4, 23))
    } else{
        plt <- plt +
            ggplot2::geom_point(alpha = 0.3, size = 1.5, shape = 4)
    }
    
    
    if(lines){
    # Plot guide lines 
        plt <- plt +
            ggplot2::geom_abline(
                slope = 1, color = "gray", linetype = "dashed") +
            ggplot2::geom_hline(
                yintercept = 0, color = "gray", linetype = "dashed") +
            ggplot2::geom_vline(
                xintercept = 0, color = "gray", linetype = "dashed")
    }
    
    # Make plot interactive if user desires
    if(plotly) plt <- plotly::ggplotly(plt, tooltip = "text")
    
    return(plt)
    }