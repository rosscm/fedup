context("Pathways")

test_that("Test that writePathways works with GMT input", {
  pathway_file <- system.file("extdata", "Human_GOBP_AllPathways_no_GO_iea_November_17_2020_symbol.gmt", package = "FEDUP")
  supported_types <- c("gmt", "txt", "xlsx")
  file_type <- sub(".*\\.", "", pathway_file)
  expect_true(file_type %in% supported_types)
  pathways <- writePathways(pathway_file, MIN_GENE = 10, MAX_GENE = 500)
  expect_true(is.list(pathways))
  expect_false(any(duplicated(names(pathways))))
})

test_that("Test that writePathways works with XLSX input", {
  pathway_file <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
  supported_types <- c("gmt", "txt", "xlsx")
  file_type <- sub(".*\\.", "", pathway_file)
  expect_true(file_type %in% supported_types)
  pathways <- writePathways(pathway_file, header = TRUE,
                            pathway_col = "Enriched.GO.names", gene_col = "Gene.ID",
                            MIN_GENE = 10, MAX_GENE = 500)
  expect_true(is.list(pathways))
  expect_false(any(duplicated(names(pathways))))
})
