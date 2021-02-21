## code to prepare `pathwaysXLSX` dataset goes here
library(tibble)

pathwayFile <- system.file("extdata", "SAFE_terms.txt", package = "fedup")
pathwaysTXT <- readPathways(
    pathwayFile,
    header = TRUE,
    pathCol = "Enriched.GO.names",
    geneCol = "Gene.ID"
)

names(pathwaysTXT) <- stringi::stri_trans_general(names(pathwaysTXT), "latin-ascii")
usethis::use_data(pathwaysTXT, compress = "xz", version = 2, overwrite = TRUE)
