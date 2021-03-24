#' Named list of human pathway annotations obtained from a GMT file.
#'
#' GMT file is available from
#' \url{http://download.baderlab.org/EM_Genesets/November_17_2020/Human/symbol}
#'
#' Raw data location
#' system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt",
#'     package = "fedup")
#' Script to prepare data
#' system.file("script", "pathwaysGMT.R", package = "fedup")
#'
#' @format a named list of 1437 vectors
#' @usage data(pathwaysGMT)
"pathwaysGMT"

#' Example list of SAFE terms obtained from a XLSX file.
#'
#' Raw data file (S5) is available from
#' \url{https://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/}
#'
#' Raw data location
#' system.file("extdata", "SAFE_terms.xlsx", package = "fedup")
#'
#' Script to prepare data
#' system.file("script", "pathwaysXLSX.R", package = "fedup")
#'
#' @format a named list of 317 vectors
#' @usage data(pathwaysXLSX)
"pathwaysXLSX"

#' Example list of SAFE terms obtained from a TXT file.
#'
#' Raw data file (S5) is available from
#' \url{https://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/}
#'
#' Raw data location
#' system.file("extdata", "SAFE_terms.txt", package = "fedup")
#'
#' Script to prepare data
#' system.file("script", "pathwaysTXT.R", package = "fedup")
#'
#' @format a named list of 317 vectors
#' @usage data(pathwaysTXT)
"pathwaysTXT"

#' Example list of a background set and single set of test genes to use for
#' enrichment analysis.
#'
#' Raw Excel data file (Sup Tab 2) is available from
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7566881/}
#'
#' Script to prepare data
#' system.file("script", "genes.R", package = "fedup")
#'
#' @format a named list with two vector elements, one common background
#' gene vector and one test gene vector.
#' @usage data(geneSingle)
"geneSingle"

#' Example list of a background set and two sets of test genes to use for
#' enrichment analysis.
#'
#' Raw Excel data file (Sup Tab 2) is available from
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7566881/}
#'
#' Script to prepare data
#' system.file("script", "genes.R", package = "fedup")
#'
#' @format a named list with three vector elements, one common background
#' gene vector and two test gene vectors.
#' @usage data(geneDouble)
"geneDouble"

#' Example list of a background set and multiple sets of test genes
#' to use for enrichment analysis.
#'
#' Raw Excel data file (Sup Tab 2) is available from
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7566881/}
#'
#' Script to prepare data
#' system.file("script", "genes.R", package = "fedup")
#'
#' @format a named list with thirteen vector elements, one common background
#' gene vector and twelve test gene vectors.
#' @usage data(geneMulti)
"geneMulti"
