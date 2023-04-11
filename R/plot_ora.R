#' Scatter plot of Over-Representation Analysis results
#' 
#' @param ora Output of \code{\link[pairedGSEA:paired_ora]{paired_ora}}
#' @param plotly (Default: \code{FALSE})
#' Logical on whether to return plot as an
#' interactive \code{\link[plotly:ggplotly]{plotly}} plot or a simple ggplot.
#' @param pattern Highlight pathways containing a specific regex pattern
#' @param cutoff (Default: \code{0.2}) Adjusted p-value cutoff for
#' pathways to include
#' @param lines (Default: \code{TRUE}) Whether to show dashed lines
#' @param colors (Default: \code{c("darkgray", "purple", "navy")})
#' Colors to use in plot. The colors are ordered as "Both", "DGS", and "DGE"
#' @importFrom stats cor
#' @importFrom S4Vectors as.data.frame
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
        = any(grepl("_splicing", colnames(ora)))
    )
    
    check_missing_package(package = "ggplot2")
    if(plotly) check_missing_package(package = "plotly")
    
    # Filter ora
    ora <- subset(
        ora, ora$padj_expression < cutoff | ora$padj_splicing < cutoff)
    stopifnot("No over-represented gene sets found." = nrow(ora) > 0)
    # Compute spearman correlation
    correlation <- spearman(ora)
    
    # Set default match to FALSE for when pattern is NULL
    ora$pattern_match <- FALSE
    if (!is.null(pattern)) {
        ora$pattern_match <- grepl(tolower(pattern), tolower(ora$pathway))
    }
    
    ora <- add_plot_color(ora, cutoff)
    # Extract matches
    matches <- ora[ora$pattern_match == TRUE,]
    
    if(!is.null(pattern) & sum(ora$pattern_match) == 0) warning(
        "No matches found with pattern: ", pattern)
    
    stopifnot(
        "Please give at least one color per category of overrepresentation" =
            length(unique(ora$plot_color)) <= length(colors))
    colors <- subset_colors(ora, colors)

    plt <- create_plot(ora)
    plt <- add_layout(plt, correlation, cutoff, colors)
    plt <- add_matches(plt, matches, pattern)
    if(lines) plt <- add_lines(plt)
    # Make plot interactive if user desires
    if(plotly) plt <- plotly::ggplotly(plt, tooltip = "text")
    
    return(plt)
}




#' Compute Spearman correlation between enrichment scores
#'
#' @inheritParams plot_ora
#'
#' @return A numeric value indicating the
#' Spearman correlation between the enrichment scores.
#' @importFrom stats cor
#' @keywords internal
#' @noRd
spearman <- function(ora) {
    correlation <- cor(
        ora$enrichment_score_splicing,
        ora$enrichment_score_expression, method = "spearman")
    return(round(correlation, 2))
}

#' Add plot color for significance
#' @noRd
add_plot_color <- function(ora, cutoff) {
    
    ora$plot_color <- ifelse(
        ora$padj_splicing < cutoff & ora$padj_expression < cutoff, "Both",
        ifelse(
            ora$padj_splicing < cutoff, "Only Splicing",
            ifelse(ora$padj_expression < cutoff, "Only Expression", "NA")))
    ora$plot_color <- factor(
        ora$plot_color, levels = c("Both", "Only Splicing", "Only Expression"))
    ora <- ora[ora$plot_color != "NA", ]
    stopifnot("No significant gene sets in input data" = nrow(ora) > 0)
    
    return(ora)
}

#' Subset colors if not all are needed
#' @noRd
subset_colors <- function(plt_data, colors) {
    return(colors[c(
        "Both", "Only Splicing", "Only Expression") %in% plt_data$plot_color])
}

#' Create aes of plt_data
#' @noRd
create_plot <- function(ora) {
    if(is(ora, "DFrame")) ora <- S4Vectors::as.data.frame(ora)
    plt <- ggplot2::ggplot(data = ora) +
        ggplot2::aes(
            x = .data$enrichment_score_splicing,
            y = .data$enrichment_score_expression,
            fill = .data$plot_color,
            color = .data$plot_color,
            text = .data$pathway)
    
    return(plt)
}

#' Add layout to plot
#' @noRd
add_layout <- function(plt, correlation, cutoff, colors){
    
    # Add correlation text to plot
    plt + 
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
        ggplot2::scale_colour_manual(values = colors)
    
}

#' Identify and add matches to plot
#' @noRd
add_matches <- function(plt, matches, pattern){
    if(is(matches, "DFrame")) matches <- S4Vectors::as.data.frame(matches)
    if(nrow(matches) > 0){
        # Add matches to plot if a pattern was given and matches found
        plt + 
            ggplot2::geom_point(alpha = 0.3, size = 1.5) +
            ggplot2::aes(
                text = paste0(.data$pathway, "\nMatch: ", .data$pattern_match),
                shape = .data$pattern_match) + 
            ggplot2::geom_point(
                data = matches, size = 2.5, alpha = 0.9, fill = "red",
                ggplot2::aes(shape = TRUE)) +
            ggplot2::labs(shape = pattern) +
            ggplot2::scale_shape_manual(values = c(4, 23))
    } else{
        plt <- plt +
            ggplot2::geom_point(alpha = 0.3, size = 1.5, shape = 4)
    }
    
    
}

#' Add guiding lines
#' @noRd
add_lines <- function(plt){
    plt +
        ggplot2::geom_abline(
            slope = 1, color = "gray", linetype = "dashed") +
        ggplot2::geom_hline(
            yintercept = 0, color = "gray", linetype = "dashed") +
        ggplot2::geom_vline(
            xintercept = 0, color = "gray", linetype = "dashed")
}
