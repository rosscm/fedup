library(biomaRt)
library(dplyr)

# Read in example pathway GMT file
pathway_file <- system.file("extdata", "Human_GOBP_AllPathways_no_GO_iea_November_17_2020_symbol.gmt", package = "FEDUP")

# Use FEDUP::writePathways() to convert to list
pathwaysGMT <- writePathways(pathway_file,
                             MIN_GENE = 10,
                             MAX_GENE = 500)

# Select olfactory signalling pathway genes to use as test genes
testGene <- pathwaysGMT[[grep("381753", names(pathwaysGMT))]]

# Grab all unique genes to use as background genes
backgroundGene <- unique(unlist(pathwaysGMT))

# Compress and save
save(testGene, file = file.path("data", "testGene.rda"), compress = "xz")
save(backgroundGene, file = file.path("data", "backgroundGene.rda"), compress = "xz")
