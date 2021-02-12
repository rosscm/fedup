context("Enrichment dotplot")

test_that("Test that plotDotPlot works", {
    data(testGene)
    data(backgroundGene)
    data(pathwaysGMT)
    fedupRes <- runFedup(testGene, backgroundGene, pathwaysGMT)
    fedupEnr <- head(fedupRes[with(fedupRes, which(status == "Enriched")), ], 10)
    fedupDep <- head(fedupRes[with(fedupRes, which(status == "Depleted")), ], 10)
    fedupPlot <- rbind(fedupEnr, fedupDep)
    fedupPlot$log10fdr <- -log10(fedupPlot$fdr + 1e-10) # log10-transform FDR for plotting
    fedupPlot$pathway <- gsub("\\%.*", "", fedupPlot$pathway) # clean pathway names
    temp <- tempfile("plot", fileext = ".png")
    png(filename = temp, width = 2750, height = 1600, res = 300)
    plotDotPlot(
        df = fedupPlot,
        xVar = "log10fdr",
        yVar = "pathway",
        xLab = "-log10(FDR)",
        fillVar = "status",
        fillLab = "Enrichment status",
        sizeVar = "fold_enrichment",
        sizeLab = "Fold enrichment"
    )
    dev.off()
    expect_true(TRUE)
})
