context("Enrichment map")

test_that("Test that writeFemap works", {
    data(testGene)
    data(backgroundGene)
    data(pathwaysGMT)
    fedupRes <- runFedup(testGene, backgroundGene, pathwaysGMT)
    resultsFile <- tempfile("fedupRes", fileext = ".txt")
    writeFemap(fedupRes, resultsFile)
    femapRes <- read.delim(resultsFile)

    expect_equal(nrow(fedupRes), nrow(femapRes))
    expect_true("status" %in% colnames(femapRes))
    expect_true(fedupRes[1, "status"] == "Enriched" && femapRes[1, "status"] == 1)
    expect_true(fedupRes[1436, "status"] == "Depleted" && femapRes[1436, "status"] == -1)
})

test_that("Test the plotFemap returns NULL when Cytoscape is not running", {
     data(testGene)
     data(backgroundGene)
     data(pathwaysGMT)
     gmtFile <- tempfile("pathwaysGMT", fileext = ".gmt")
     fedupRes <- runFedup(testGene, backgroundGene, pathwaysGMT)
     resultsFile <- tempfile("fedupRes", fileext = ".txt")
     netFile <- tempfile("fedup_EM", fileext = ".png")
     writePathways(pathwaysGMT, gmtFile)
     writeFemap(fedupRes, resultsFile)
     expect_null(plotFemap(
         gmtFile = gmtFile,
         resultsFile = resultsFile,
         qvalue = 0.05,
         netName = "fedup_EM",
         netFile = netFile
     ))
})
