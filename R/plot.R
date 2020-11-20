#' Visualizes pathway enrichment and depletion using ggplot.
#' @param df (data.frame) table with FEDUP enrichment results to plot
#' @param x_var (char) x-axis variable (must be a column value in `df`)
#' @param y_var (char) y-axis variable (must be a column value in `df`)
#' @param x_lab (char) x-axis label (default `x_var`)
#' @param y_lab (char) y-axis label (default NULL)
#' @param p_title (char) plot title (default NULL)
#' @param fill_var (char) point fill variable (default NULL)
#' @param fill_col (char) point fill colours (default NULL)
#' @param fill_lab (char) point fill label (default `fill_var`)
#' @param size_var (char) point size variable (default NULL)
#' @param size_lab (char) point size label (default `size_var`)
#' @return ggplot object with the enrichment dot plot
#' @import ggplot2
#' @import ggthemes
#' @import forcats
#' @import RColorBrewer
#' @export
#' @examples
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' fedup_res$log10fdr <- -log10(fedup_res$fdr) # log10-transform FDR
#' fedup_plot <- fedup_res[head(order(fdr, real_pathway_frac), 20),]
#' \dontrun{
#' plotDotPlot(df = fedup_plot,
#'             x_var = "log10fdr",
#'             y_var = "pathways",
#'             x_lab = "-log10(FDR)",
#'             fill_var = "enrichment",
#'             fill_lab = "Enrichment",
#'             size_var = "real_pathway_frac",
#'             size_lab = "Gene fraction")
#' }
plotDotPlot <- function(df,
                        x_var,
                        y_var,
                        x_lab = x_var,
                        y_lab = NULL,
                        p_title = NULL,
                        fill_var = NULL,
                        fill_col = NULL,
                        fill_lab = fill_var,
                        size_var = NULL,
                        size_lab = size_var) {

  # Set point fill if `fill_var` is specified but `fill_col` is not
  if (!is.null(fill_var) && is.null(fill_col)) {
    fill_n <- length(unique(df[[fill_var]]))
    pal_n <- ifelse(fill_n >= 3, fill_n, 3)
    fill_col <- brewer.pal(pal_n, "Pastel2")
  }

  # Plot dot plot with specified parameters
  p <- ggplot(df, aes_string(x = x_var,
                             y = fct_reorder(df[[y_var]], df[[x_var]]),
                             fill = fill_var,
                             size = size_var)) +
          geom_point(shape = 21, colour = "black") +
          labs(x = x_lab,
               y = y_lab,
               title = p_title,
               fill = fill_lab,
               size = size_lab) +
          scale_fill_manual(values = fill_col) +
          theme_few(base_size = 14) +
          theme(plot.title = element_text(hjust = 0.5),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.key.size = unit(0.1, "line"),
                panel.grid.major.y = element_line(linetype = "dotted", colour = "lightgrey"),
                aspect.ratio = 16/9)

  if (is.numeric(df[[x_var]])) {
    # Set x-axis limits so points are not cut off from plot window
    xmin <- floor(min(df[[x_var]]))
    xmax <- ceiling(max(df[[x_var]]))
    p <- p + xlim(xmin, xmax)
  }
  return(p)
}
