library(openxlsx)
library(tibble)

# Input XLSX file
pathway_file <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")

# Use FEDUP::writePathways() to convert to list
pathwaysXLSX <- writePathways(pathway_file,
                              header = TRUE,
                              pathway_col = "Enriched.GO.names",
                              gene_col = "Gene.ID",
                              MIN_GENE = 10,
                              MAX_GENE = 500)

# Compress and save
save(pathwaysXLSX, file = file.path("..", "data", "pathwaysXLSX.rda"), compress = "xz")
