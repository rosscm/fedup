## code to prepare `testGene` and `backgroundGene` datasets goes here
pathwayFile <- system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt", package = "FEDUP")
pathwaysGMT <- readPathways(pathwayFile, minGene = 10, maxGene = 500)

testGene <- pathwaysGMT[[grep("397014", names(pathwaysGMT))]] # Reactome muscle contraction pathway
backgroundGene <- unique(unlist(pathwaysGMT))

usethis::use_data(testGene, compress = "xz", version = 2, overwrite = TRUE)
usethis::use_data(backgroundGene, compress = "xz", version = 2, overwrite = TRUE)
