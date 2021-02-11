## code to prepare `pathwaysXLSX` dataset goes here
library(tibble)

pathway_file <- system.file("extdata", "SAFE_terms.txt", package = "FEDUP")
pathwaysTXT <- readPathways(
  pathway_file,
  header = TRUE,
  pathway_col = "Enriched.GO.names",
  gene_col = "Gene.ID")

names(pathwaysTXT) <- stringi::stri_trans_general(names(pathwaysTXT), "latin-ascii")
usethis::use_data(pathwaysTXT, compress = "xz", version = 2, overwrite = TRUE)
