## code to prepare `pathwaysXLSX` dataset goes here
library(openxlsx)
library(tibble)

pathwayFile <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
pathwaysXLSX <- readPathways(
    pathwayFile,
    header = TRUE,
    pathCol = "Enriched.GO.names",
    geneCol = "Gene.ID"
)

names(pathwaysXLSX) <- stringi::stri_trans_general(names(pathwaysXLSX), "latin-ascii")
usethis::use_data(pathwaysXLSX, compress = "xz", version = 2, overwrite = TRUE)
