# Download raw data from https://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/
# Data file S5 (sheet 3)
library(openxlsx)
library(tibble)
library(biomaRt)
library(dplyr)

# Use biomaRt to get gene members per SAFE term
#pathwayFile <- system.file("extdata", "Data_File_S5_SAFE_analysis_Gene_cluster_identity_and_functional_enrichments.xlsx", package = "fedup")
#pathway <- read.xlsx(pathwayFile, sheet = 3)

# Query Ensembl for gene symbols annotated to SAFE terms
#ensembl <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
#ensembl_gene <- getBM(
#    attributes = c("go_id", "ensembl_gene_id", "external_gene_name"),
#    mart = ensembl
#)
#colnames(ensembl_gene) <- c("Enriched.GO.IDs", "ORF.ID", "Gene.ID")
#pathway <- left_join(pathway, ensembl_gene, by = "Enriched.GO.IDs")
#write.xlsx(pathway, file.path("inst", "extdata", "SAFE_terms.xlsx"), row.names = FALSE)

# Raw data file annotated with gene symbols
pathwayFile <- system.file("extdata", "SAFE_terms.xlsx", package = "fedup")
pathwaysXLSX <- readPathways(
    pathwayFile,
    header = TRUE,
    pathCol = "Enriched.GO.names",
    geneCol = "Gene.ID"
)

names(pathwaysXLSX) <- stringi::stri_trans_general(names(pathwaysXLSX), "latin-ascii")
usethis::use_data(pathwaysXLSX, compress = "xz", version = 2, overwrite = TRUE)
