## code to prepare `pathwaysGMT` dataset goes here
pathwayFile <- system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt", package = "fedup")
pathwaysGMT <- readPathways(pathwayFile, minGene = 10, maxGene = 500)
names(pathwaysGMT) <- stringi::stri_trans_general(names(pathwaysGMT), "latin-ascii")
usethis::use_data(pathwaysGMT, compress = "xz", version = 2, overwrite = TRUE)
