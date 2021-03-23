#' Prepares input gene list for \link[fedup]{runFedup}.
#'
#' @description This function takes any number of test genes and a common
#' background set of genes and properly formats them for to pass to
#' \link[fedup]{runFedup} \code{gene} argument.
#' @param setName (char) character vector naming gene vectors passed to
#' \code{...} (must include "background" (case sensitive)).
#' @param ... (char) \code{n} vectors with at least one background set of genes
#' and any number of test genes.
#' @return List of length \code{n} with gene vectors corresponding to those
#' passed to \code{...}.
#' @examples
#' # Raw gene data file
#' genesFile <- system.file("extdata",
#'     "NIHMS1587165-supplement-1587165_Sup_Tab_2.txt", package = "fedup")
#' genes <- read.delim(genesFile, h = TRUE, as.is = TRUE)
#' # Prepare gene vectors
#' b <- unique(genes[, "gene"])
#' set1 <- unique(genes[which(genes$FASN_merge < -0.4), "gene"])
#' set2 <- unique(genes[which(genes$FASN_merge > 0.4), "gene"])
#' setName <- c("background", "negative", "positive")
#' geneDouble <- prepInput(setName, b, set1, set2)
#' @export
prepInput <- function(setName, ...) {
    inList <- list(...)
    if (length(inList) < 2) {
        stop("Input gene list must include at least one background and one test
        set of genes to run enrichment analysis.")
    }
    if (FALSE %in% vapply(inList, is.vector, logical(1))) {
        stop("Objects passed through '...' must be in vector format.")
    }
    if (length(setName) != length(inList)) {
        stop("Length of 'setName' does not match with number of inputs.
        Make sure all inputs are named via 'setName'.")
    }
    if (length(which(duplicated(setName)))) {
        stop("Names provided to the 'setName' vector must be unique.")
    }
    if (!"background" %in% setName) {
        stop("Input gene list must have a named background set. Make sure the
        'setName' vector includes 'background' (case sensitive!).")
    }
    names(inList) <- setName
    return(inList)
}
