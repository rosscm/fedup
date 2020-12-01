# Input GMT file
# Downloaded from http://download.baderlab.org/EM_Genesets/November_17_2020/Human/symbol/Human_GOBP_AllPathways_no_GO_iea_November_17_2020_symbol.gmt
pathway_file <- system.file("extdata", "Human_GOBP_AllPathways_no_GO_iea_November_17_2020_symbol.gmt", package = "FEDUP")

# Use FEDUP::writePathways() to convert to list
pathwaysGMT <- writePathways(pathway_file,
                             MIN_GENE = 10,
                             MAX_GENE = 500)

# Compress and save
save(pathwaysGMT, file = file.path("data", "pathwaysGMT.rda"), compress = "xz")
