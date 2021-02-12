#' Visualizes pathway enrichment and depletion using ggplot.
#'
#' @param df (data.frame) table with FEDUP enrichment results to plot.
#' @param xVar (char) x-axis variable (must be a column value in \code{df}).
#' @param yVar (char) y-axis variable (must be a column value in \code{df}).
#' @param xLab (char) x-axis label (default \code{xVar} value).
#' @param yLab (char) y-axis label (default NULL).
#' @param pTitle (char) plot title (default NULL).
#' @param fillVar (char) point fill variable (default NULL).
#' @param fillCol (char) point fill colours (default NULL).
#' @param fillLab (char) point fill label (default \code{fillVar} value).
#' @param sizeVar (char) point size variable (default NULL).
#' @param sizeLab (char) point size label (default \code{sizeVar} value).
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
#'     df=fedup_plot,
#'     xVar="log10qvalue",
#'     yVar="pathway",
#'     xLab="-log10(Qvalue)",
#'     fillVar="status",
#'     fillLab="Enrichment\nstatus",
#'     sizeVar="fold_enrichment",
#'     sizeLab="Fold enrichment")
#' @import ggplot2
#' @importFrom ggthemes theme_clean
#' @importFrom forcats fct_reorder
#' @importFrom RColorBrewer brewer.pal
#' @export
plotDotPlot <- function(df, xVar, yVar,
                        xLab=xVar, yLab=NULL, pTitle=NULL,
                        fillVar=NULL, fillCol=NULL, fillLab=fillVar,
                        sizeVar=NULL, sizeLab=sizeVar) {

    if (!is.null(fillVar) && is.null(fillCol)) {
        fill_n <- length(unique(df[[fillVar]]))
        pal_n <- ifelse(fill_n >= 3, fill_n, 3)
        fillCol <- brewer.pal(pal_n, "Set1")
    }

    if (fillVar == "status") {
        df[[fillVar]] <- factor(
            df[[fillVar]],
            levels=c("Enriched", "Depleted"))
    }

    p <- ggplot(df, aes_string(
                x=xVar,
                y=fct_reorder(df[[yVar]], df[[xVar]]),
                fill=fillVar,
                size=sizeVar)) +
        geom_point(shape=21, colour="black") +
        labs(x=xLab, y=yLab, title=pTitle,
            fill=fillLab, size=sizeLab) +
        scale_fill_manual(values=fillCol) +
        theme_clean(base_size=10) +
        theme(plot.title=element_text(hjust=0.5),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10),
            legend.key.size=unit(0.1, "line"),
            plot.background=element_blank())

    # Increase x-axis limits to keep points in plot window
    if (is.numeric(df[[xVar]])) {
        xmin <- floor(min(df[[xVar]]))
        xmax <- ceiling(max(df[[xVar]]))
        p <- p + xlim(xmin, xmax)
    }
    return(p)
}
