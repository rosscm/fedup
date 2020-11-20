# Input GMT file
# Downloaded from http://ge-lab.org/gskb/3-yeast-resiluts-database/YeastDatabase_GO_gmt.gmt
pathway_file <- system.file("extdata", "YeastDatabase_GO_gmt.gmt", package = "FEDUP")

# Use FEDUP::writePathways() to convert to list
pathwaysGMT <- writePathways(pathway_file,
                             MIN_GENE = 10,
                             MAX_GENE = 500,
                             GO_class = "BP")

# Compress and save
save(pathwaysGMT, file = file.path("..", "data", "pathwaysGMT.rda"), compress = "xz")
