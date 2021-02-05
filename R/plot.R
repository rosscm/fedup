#' Visualizes pathway enrichment and depletion using ggplot.
#'
#' @param df (data.frame) table with FEDUP enrichment results to plot.
#' @param x_var (char) x-axis variable (must be a column value in `df`).
#' @param y_var (char) y-axis variable (must be a column value in `df`).
#' @param x_lab (char) x-axis label (default `x_var` value).
#' @param y_lab (char) y-axis label (default NULL).
#' @param p_title (char) plot title (default NULL).
#' @param fill_var (char) point fill variable (default NULL).
#' @param fill_col (char) point fill colours (default NULL).
#' @param fill_lab (char) point fill label (default `fill_var` value).
#' @param size_var (char) point size variable (default NULL).
#' @param size_lab (char) point size label (default `size_var` value).
#' @return object returned from ggplot with the enrichment dot plot.
#' @examples
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' fedup_plot <- fedup_res[which(fedup_res$qvalue < 0.05),]
#' fedup_plot$log10qvalue <- -log10(fedup_plot$qvalue + 1e-10)
#' fedup_plot$pathway <- gsub("\\%.*", "", fedup_plot$pathway)
#' plotDotPlot(
#'      df = fedup_plot,
#'      x_var = "log10qvalue",
#'      y_var = "pathway",
#'      x_lab = "-log10(Qvalue)",
#'      fill_var = "enrichment",
#'      fill_lab = "Enrichment\nstatus",
#'      size_var = "real_pathway_frac",
#'      size_lab = "Gene fraction")
#' @import ggplot2
#' @importFrom ggthemes theme_clean
#' @importFrom forcats fct_reorder
#' @importFrom RColorBrewer brewer.pal
#' @export
plotDotPlot <- function(df, x_var, y_var,
                        x_lab = x_var, y_lab = NULL, p_title = NULL,
                        fill_var = NULL, fill_col = NULL, fill_lab = fill_var,
                        size_var = NULL, size_lab = size_var) {
    # Set point fill if `fill_var` is specified but `fill_col` is not
    if (!is.null(fill_var) && is.null(fill_col)) {
        fill_n <- length(unique(df[[fill_var]]))
        pal_n <- ifelse(fill_n >= 3, fill_n, 3)
        fill_col <- brewer.pal(pal_n, "Set1")
    }

    # Set factor level if `fill_var` is set to `enrichment`
    if (fill_var == "enrichment") {
        df[[fill_var]] <- factor(
            df[[fill_var]],
            levels = c("Enriched", "Depleted")
        )
    }

    # Plot dot plot with specified parameters
    p <- ggplot(df, aes_string(
                x = x_var,
                y = fct_reorder(df[[y_var]], df[[x_var]]),
                fill = fill_var,
                size = size_var)) +
        geom_point(shape = 21, colour = "black") +
        labs(x = x_lab, y = y_lab, title = p_title,
            fill = fill_lab, size = size_lab) +
        scale_fill_manual(values = fill_col) +
        theme_clean(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.1, "line"),
            plot.background = element_blank())

    if (is.numeric(df[[x_var]])) {
        # Set x-axis limits so points are not cut off from plot window
        xmin <- floor(min(df[[x_var]]))
        xmax <- ceiling(max(df[[x_var]]))
        p <- p + xlim(xmin, xmax)
    }
    return(p)
}
