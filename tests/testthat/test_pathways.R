context("Pathways")

test_that("Test that readPathways stops without proper inputs", {
    expect_error(readPathways("test.123.xls"))
    expect_error(readPathways("test.gmt.123"))

    pathwayFile <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
    expect_error(readPathways(
        pathwayFile,
        header = TRUE,
        pathCol = "Enriched.GO.names", geneCol = "oops"
    ))
    expect_error(readPathways(
        pathwayFile,
        header = TRUE,
        pathCol = "oops", geneCol = "Gene.ID"
    ))
    expect_error(readPathways(
        pathwayFile,
        header = TRUE, minGene = 500,
        pathCol = "Enriched.GO.names", geneCol = "Gene.ID"
    ))
})

test_that("Test that readPathways works with GMT input", {
    pathwayFile <- system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt", package = "FEDUP")
    s <- c("gmt", "txt", "xlsx")
    f <- sub(".*\\.", "", pathwayFile)
    expect_true(f %in% s)

    pathways <- readPathways(pathwayFile, minGene = 10, maxGene = 500)
    expect_true(is.list(pathways))
    expect_equal(length(pathways), 1437)
    expect_equal(length(readPathways(pathwayFile, minGene = 10, maxGene = 500, header = TRUE)), 1436)
})

test_that("Test that readPathways works with XLSX input", {
    pathwayFile <- system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP")
    s <- c("gmt", "txt", "xlsx")
    f <- sub(".*\\.", "", pathwayFile)
    expect_true(f %in% s)

    pathways <- readPathways(
        pathwayFile,
        header = TRUE,
        pathCol = "Enriched.GO.names", geneCol = "Gene.ID"
    )
    expect_true(is.list(pathways))
    expect_equal(length(pathways), 30)
})

test_that("Test that readPathways works with TXT input", {
    pathwayFile <- system.file("extdata", "SAFE_terms.txt", package = "FEDUP")
    s <- c("gmt", "txt", "xlsx")
    f <- sub(".*\\.", "", pathwayFile)
    expect_true(f %in% s)

    pathways <- readPathways(
        pathwayFile,
        header = TRUE,
        pathCol = "Enriched.GO.names", geneCol = "Gene.ID"
    )
    expect_true(is.list(pathways))
    expect_equal(length(pathways), 30)
})

test_that("Test that writePathways works", {
    data(pathwaysXLSX)
    gmtFile <- tempfile("pathwaysXLSX", fileext = ".gmt")

    writePathways(pathwaysXLSX, gmtFile)
    pathways <- readPathways(gmtFile, header = FALSE)

    expect_equal(length(pathwaysXLSX), length(pathways))
    expect_true(is.list(pathways))
})
