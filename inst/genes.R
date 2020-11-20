library(biomaRt)
library(dplyr)

# Query S.cerevisiae GO terms and associated gene symbols
ensembl <- useMart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
go_gene <- getBM(attributes = c("go_id", "external_gene_name"), mart = ensembl)
go_gene[go_gene == ""] <- NA
go_gene <- na.omit(go_gene)

# Select random GO term and grab all genes in term to use as test genes
set.seed(10)
go_rand <- sample(go_gene$go_id, 1)
testGene <- go_gene %>%
  filter(go_id == go_rand) %>%
  select(external_gene_name) %>%
  na.omit %>%
  unlist %>%
  as.character

# Grab all unique genes to use as background genes
backgroundGene <- as.character(unique(go_gene$external_gene_name))

# Compress and save
save(testGene, file = file.path("..", "data", "testGene.rda"), compress = "xz")
save(backgroundGene, file = file.path("..", "data", "backgroundGene.rda"), compress = "xz")
