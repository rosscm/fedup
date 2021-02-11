context("Pathways")

test_that("Test that readPathways stops without proper inputs", {
    expect_error(readPathways("test.123.xls"))
    expect_error(readPathways("test.gmt.123"))

    pathway_file <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
    expect_error(readPathways(
        pathway_file, header = TRUE,
        pathway_col = "Enriched.GO.names", gene_col = "oops"))
    expect_error(readPathways(
        pathway_file, header = TRUE,
        pathway_col = "oops", gene_col = "Gene.ID"))
    expect_error(readPathways(
        pathway_file, header = TRUE, MIN_GENE = 500,
        pathway_col = "Enriched.GO.names", gene_col = "Gene.ID"))
})

test_that("Test that readPathways works with GMT input", {
    pathway_file <- system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt", package = "FEDUP")
    supported_types <- c("gmt", "txt", "xlsx")
    file_type <- sub(".*\\.", "", pathway_file)
    expect_true(file_type %in% supported_types)

    pathways <- readPathways(pathway_file, MIN_GENE = 10, MAX_GENE = 500)
    expect_true(is.list(pathways))
    expect_equal(length(pathways), 1437)
    expect_false(any(duplicated(names(pathways))))
    expect_equal(length(readPathways(pathway_file, MIN_GENE = 10, MAX_GENE = 500, header = TRUE)), 1436)
})

test_that("Test that readPathways works with XLSX input", {
    pathway_file <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
    supported_types <- c("gmt", "txt", "xlsx")
    file_type <- sub(".*\\.", "", pathway_file)
    expect_true(file_type %in% supported_types)

    pathways <- readPathways(
        pathway_file, header = TRUE,
        pathway_col = "Enriched.GO.names", gene_col = "Gene.ID")
    expect_true(is.list(pathways))
    expect_equal(length(pathways), 30)
    expect_false(any(duplicated(names(pathways))))
})

test_that("Test that readPathways works with TXT input", {
    pathway_file <- system.file("extdata", "SAFE_terms.txt", package = "FEDUP")
    supported_types <- c("gmt", "txt", "xlsx")
    file_type <- sub(".*\\.", "", pathway_file)
    expect_true(file_type %in% supported_types)

    pathways <- readPathways(
        pathway_file, header = TRUE,
        pathway_col = "Enriched.GO.names", gene_col = "Gene.ID")
    expect_true(is.list(pathways))
    expect_equal(length(pathways), 30)
    expect_false(any(duplicated(names(pathways))))
})

test_that("Test that writePathways works", {
    data(pathwaysXLSX)
    gmt_file <-  tempfile("pathwaysXLSX", fileext = ".gmt")

    writePathways(pathwaysXLSX, gmt_file)
    pathways <- readPathways(gmt_file, header = FALSE)

    expect_equal(length(pathwaysXLSX), length(pathways))
    expect_true(is.list(pathways))
})
