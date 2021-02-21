#' Named list of human pathway annotations obtained from a GMT file.
#'
#' Raw GMT file is available from
#' \url{http://download.baderlab.org/EM_Genesets/November_17_2020/Human/symbol}
#'
#' Raw data location
#' system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt",
#'     package = "fedup")
#' Script to prepare data
#' system.file("data-raw", "pathwaysGMT.R", package ="fedup")
#'
#' @format a named list of 1437 vectors
"pathwaysGMT"

#' Example list of yeast SAFE terms obtained from a XLSX file.
#'
#' Raw data location
#' system.file("extdata", "SAFE_terms.xlsx", package = "fedup")
#'
#' Script to prepare data
#' system.file("data-raw", "pathwaysXLSX.R", package = "fedup")
#'
#' @format a named list of 30 vectors
"pathwaysXLSX"

#' Example list of yeast SAFE terms obtained from a TXT file.
#'
#' Raw data location
#' system.file("extdata", "SAFE_terms.txt", package = "fedup")
#'
#' Script to prepare data
#' system.file("data-raw", "pathwaysTXT.R", package = "fedup")
#'
#' @format a named list of 30 vectors
"pathwaysTXT"

#' Example vector of human genes to use as test set for enrichment.
#'
#' Script to prepare data
#' system.file("data-raw", "genes.R", package = "fedup")
#'
#' @format a character vector with 190 elements (gene IDs)
"testGene"

#' Example vector of human genes to use as background set for enrichment.
#'
#' Script to generate data
#' system.file("data-raw", "genes.R", package = "fedup")
#'
#' @format a character vector with 10208 elements (gene IDs)
"backgroundGene"
