#' Visualizes fedup enrichment and depletion results using ggplot.
#'
#' @description This function supports any combination of numeric x-y variables
#' to plot from fedup results. The list outputted by \link[fedup]{runFedup}
#' must first be converted to a data.frame before plotting (see examples for
#' sample use).
#' @param df (data.frame) table with fedup results generated via
#' \link[fedup]{runFedup}
#' @param xVar (char) x-axis variable (must be a column value in \code{df})
#' @param yVar (char) y-axis variable (must be a column value in \code{df})
#' @param xLab (char) x-axis label (default \code{xVar} value)
#' @param yLab (char) y-axis label (default NULL)
#' @param pTitle (char) plot title (default NULL)
#' @param fillVar (char) point fill variable (default NULL)
#' @param fillCol (char) point fill colours (default NULL)
#' @param fillLab (char) point fill label (default \code{fillVar} value)
#' @param sizeVar (char) point size variable (default NULL)
#' @param sizeLab (char) point size label (default \code{sizeVar} value)
#' @return Object returned from ggplot with the enrichment dot plot.
#' @examples
#' # Load example data
#' data(geneDouble)
#' data(pathwaysGMT)
#' # Load external libraries
#' suppressMessages(library(dplyr))
#' suppressMessages(library(tidyr))
#' # Run fedup
#' fedupRes <- runFedup(geneDouble, pathwaysGMT)
#' # Prepare dataframe from fedup results
#' fedupPlot <- fedupRes %>%
#'     bind_rows(.id = "set") %>%
#'     separate(col = "set", into = c("set", "sign"), sep = "_") %>%
#'     subset(qvalue < 0.01) %>%
#'     mutate(log10qvalue = -log10(qvalue)) %>%
#'     mutate(pathway = gsub("\\%.*", "", pathway)) %>%
#'     as.data.frame()
#' # Plot
#' p <- plotDotPlot(
#'     df = fedupPlot,
#'     xVar = "log10qvalue",
#'     yVar = "pathway",
#'     xLab = "-log10(qvalue)",
#'     fillVar = "sign",
#'     fillLab = "Genetic interaction",
#'     fillCol = c("#0077f1", "#fcde24"),
#'     sizeVar = "fold_enrichment",
#'     sizeLab = "Fold enrichment"
#' )
#' @import ggplot2
#' @importFrom ggthemes theme_clean
#' @importFrom forcats fct_reorder
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @usage
#' plotDotPlot(
#'     df,
#'     xVar,
#'     yVar,
#'     xLab = xVar,
#'     yLab = NULL,
#'     pTitle = NULL,
#'     fillVar = NULL,
#'     fillCol = NULL,
#'     fillLab = fillVar,
#'     sizeVar = NULL,
#'     sizeLab = sizeVar
#' )
plotDotPlot <- function(df, xVar, yVar,
    xLab = xVar, yLab = NULL, pTitle = NULL,
    fillVar = NULL, fillCol = NULL, fillLab = fillVar,
    sizeVar = NULL, sizeLab = sizeVar) {

    if (!is.null(fillVar) && is.null(fillCol)) {
        fill_n <- length(unique(df[[fillVar]]))
        pal_n <- ifelse(fill_n >= 3, fill_n, 3)
        fillCol <- brewer.pal(pal_n, "Set1")
    }

    if (is.numeric(df[[xVar]])) {
        p <- ggplot(df, aes_string(
            x = xVar,
            y = fct_reorder(df[[yVar]], df[[xVar]])
        ))
    } else {
        p <- ggplot(df, aes_string(x = xVar, y = yVar))
    }
    p <- p +
        geom_point(
            aes_string(fill = fillVar, size = sizeVar),
            shape = 21,
            colour = "black"
        ) +
        labs(
            x = xLab, y = yLab, title = pTitle,
            fill = fillLab, size = sizeLab
        ) +
        scale_fill_manual(values = fillCol) +
        theme_clean(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.1, "line"),
            plot.background = element_blank()
        )
    return(p)
}
