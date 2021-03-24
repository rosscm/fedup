# Download raw data from
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7566881/bin/NIHMS1587165-supplement-1587165_Sup_Tab_2.xlsx
# Supplementary Table 2 (sheet 1)
library(openxlsx)

# Reformat data file
genesFile <- system.file("extdata", "NIHMS1587165-supplement-1587165_Sup_Tab_2.txt", package = "fedup")
genes <- read.delim(genesFile, h = TRUE, as.is = TRUE)

# Prepare test and background gene data
## geneSingle (one test set vs. background)
geneSingle <- list()
geneSingle[["background"]] <- unique(genes[,"gene"])
geneSingle[["FASN_negative"]] <- unique(genes[which(genes$FASN_merge < -0.4), "gene"])
usethis::use_data(geneSingle, compress = "xz", version = 2, overwrite = TRUE)

## geneDouble (two test sets vs. background)
geneDouble <- list()
geneDouble[["background"]] <- unique(genes[,"gene"])
geneDouble[["FASN_negative"]] <- unique(genes[which(genes$FASN_merge < -0.4), "gene"])
geneDouble[["FASN_positive"]] <- unique(genes[which(genes$FASN_merge > 0.4), "gene"])
usethis::use_data(geneDouble, compress = "xz", version = 2, overwrite = TRUE)

## geneMulti (multiple test sets vs. background)
geneMulti <- list()
geneMulti[["background"]] <- unique(genes[,"gene"])
geneMulti[["FASN_negative"]] <- unique(genes[which(genes$FASN_merge < -0.4), "gene"])
geneMulti[["FASN_positive"]] <- unique(genes[which(genes$FASN_merge > 0.4), "gene"])
geneMulti[["ACACA_negative"]] <- unique(genes[which(genes$ACACA < -0.4), "gene"])
geneMulti[["ACACA_positive"]] <- unique(genes[which(genes$ACACA > 0.4), "gene"])
geneMulti[["LDLR_negative"]] <- unique(genes[which(genes$LDLR < -0.4), "gene"])
geneMulti[["LDLR_positive"]] <- unique(genes[which(genes$LDLR > 0.4), "gene"])
geneMulti[["SREBF1_negative"]] <- unique(genes[which(genes$SREBF1 < -0.4), "gene"])
geneMulti[["SREBF1_positive"]] <- unique(genes[which(genes$SREBF1 > 0.4), "gene"])
geneMulti[["SREBF2_negative"]] <- unique(genes[which(genes$SREBF2 < -0.4), "gene"])
geneMulti[["SREBF2_positive"]] <- unique(genes[which(genes$SREBF2 > 0.4), "gene"])
geneMulti[["C12orf49_negative"]] <- unique(genes[which(genes$C12orf49 < -0.4), "gene"])
geneMulti[["C12orf49_positive"]] <- unique(genes[which(genes$C12orf49 > 0.4), "gene"])
usethis::use_data(geneMulti, compress = "xz", version = 2, overwrite = TRUE)
