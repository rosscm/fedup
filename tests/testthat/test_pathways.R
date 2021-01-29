context("Pathways")

test_that("Test that readPathways works with GMT input", {
  pathway_file <- system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt", package = "FEDUP")
  supported_types <- c("gmt", "txt", "xlsx")
  file_type <- sub(".*\\.", "", pathway_file)
  expect_true(file_type %in% supported_types)
  pathways <- readPathways(pathway_file, MIN_GENE = 10, MAX_GENE = 500)
  expect_true(is.list(pathways))
  expect_false(any(duplicated(names(pathways))))
})

test_that("Test that readPathways works with XLSX input", {
  pathway_file <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
  supported_types <- c("gmt", "txt", "xlsx")
  file_type <- sub(".*\\.", "", pathway_file)
  expect_true(file_type %in% supported_types)
  pathways <- readPathways(pathway_file, header = TRUE,
                           pathway_col = "Enriched.GO.names", gene_col = "Gene.ID",
                           MIN_GENE = 10, MAX_GENE = 500)
  expect_true(is.list(pathways))
  expect_false(any(duplicated(names(pathways))))
})
