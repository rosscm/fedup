## code to prepare `pathwaysGMT` dataset goes here
pathway_file <- system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt", package = "FEDUP")
pathwaysGMT <- readPathways(pathway_file, MIN_GENE = 10, MAX_GENE = 500)
names(pathwaysGMT) <- stringi::stri_trans_general(names(pathwaysGMT), "latin-ascii")
usethis::use_data(pathwaysGMT, compress = "xz", version = 2, overwrite = TRUE)
