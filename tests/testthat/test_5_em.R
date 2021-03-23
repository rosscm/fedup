context("(5) Network")

test_that("Test that writeFemap works", {
    data(geneDouble)
    data(pathwaysGMT)
    fedupRes <- runFedup(geneDouble, pathwaysGMT)
    resultsFolder <- tempdir()
    writeFemap(fedupRes, resultsFolder)
    femapFiles <- list.files(pattern = "femap.*.txt", path = resultsFolder, full.names = TRUE)
    femapRes <- lapply(femapFiles, read.delim)

    expect_equal(length(fedupRes), length(femapRes))
    expect_equal(nrow(fedupRes)[[1]], nrow(fedupRes)[[1]])
    expect_true("status" %in% colnames(femapRes[[1]]))
    expect_true(fedupRes[[1]][1, "status"] == "enriched" && femapRes[[1]][1, "status"] == 1)
    expect_true(fedupRes[[1]][1437, "status"] == "depleted" && femapRes[[1]][1436, "status"] == -1)
})

test_that("Test the plotFemap returns NULL when Cytoscape is not running", {
    data(geneDouble)
    data(pathwaysGMT)

    fedupRes <- runFedup(geneDouble, pathwaysGMT)
    resultsFolder <- tempdir()
    writeFemap(fedupRes, resultsFolder)
    gmtFile <- tempfile("pathwaysGMT", fileext = ".gmt")
    writePathways(pathwaysGMT, gmtFile)
    netFile <- tempfile("fedup_EM", fileext = ".png")
    expect_null(plotFemap(
        gmtFile = gmtFile,
        resultsFolder = resultsFolder,
        qvalue = 0.05,
        netName = "fedup_EM",
        netFile = netFile
    ))
})
