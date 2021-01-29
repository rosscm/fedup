## code to prepare `pathwaysXLSX` dataset goes here
library(openxlsx)
library(tibble)

pathway_file <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
pathwaysXLSX <- readPathways(
  pathway_file,
  header = TRUE,
  pathway_col = "Enriched.GO.names",
  gene_col = "Gene.ID",
  MIN_GENE = 10,
  MAX_GENE = 500
)

names(pathwaysXLSX) <- stringi::stri_trans_general(names(pathwaysXLSX), "latin-ascii")
usethis::use_data(pathwaysXLSX, compress = "xz", version = 2, overwrite = TRUE)
