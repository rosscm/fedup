#' # Prepares fedup input.
#' @param test_genes (char)
#' @param background_genes (char)
#' @param pathways (char)
#' @export
inputObject <- function(test_genes, background_genes, pathways) {
  # Stop if test_genes is empty

}


#' Runs gene enrichment and depletion analysis for a list of pathways.
#' @param test_genes (char)
#' @param background_genes (char)
#' @param pathways (char) list of gene sets to test for enrichment
#' @return table of gene enrichment and depletion results. Rows represent tested
#' pathways. Columns represent:
#' \itemize{
#'    \item pathway --
#'    \item
#'    \item
#'}
#' @importFrom utils head read.delim tail
#' @export
runFedup <- function(test_genes, background_genes, pathways) {

}
