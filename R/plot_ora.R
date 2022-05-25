


plot_ora <- function(ora, plotly = TRUE){
  
  check_missing_package(package = "ggplot2")
  if(plotly) check_missing_package(package = "plotly")
  
  plt <- ora %>% 
    filter(padj_deseq < 0.05 | padj_dexseq < 0.05) %>% 
    ggplot2::ggplot() +
    ggplot2::aes(x = log2(relative_risk_dexseq+0.06),
                 y = log2(relative_risk_deseq+0.06),
                 text = pathway) +
    ggplot2::geom_point() + 
    ggplot2::geom_abline(slope = 1, color = "gray", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = c(-4, 6), ylim = c(-4,6)) +
    ggplot2::theme_classic()
  if(plotly) plt <- plotly::ggplotly(plt)
  
  return(plt)
}