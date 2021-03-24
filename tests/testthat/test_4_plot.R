context("(4) Plot")

data(geneDouble)
data(pathwaysGMT)
fedupRes <- runFedup(geneDouble, pathwaysGMT)
fedupPlot <- fedupRes %>%
    dplyr::bind_rows(.id = "set") %>%
    tidyr::separate(col = "set", into = c("set", "sign"), sep = "_") %>%
    subset(qvalue < 0.01) %>%
    mutate(log10qvalue = -log10(qvalue)) %>%
    mutate(pathway = gsub("\\%.*", "", pathway)) %>%
    as.data.frame()

test_that("Test that plotDotPlot works with numeric class x variable", {
    temp <- tempfile("plot", fileext = ".png")
    png(filename = temp, width = 2750, height = 1600, res = 300)
    plotDotPlot(
        df = fedupPlot,
        xVar = "log10qvalue",
        yVar = "pathway",
        xLab = "-log10(qvalue)",
        fillVar = "sign",
        fillLab = "Genetic interaction",
        sizeVar = "fold_enrichment",
        sizeLab = "Fold enrichment"
    )
    dev.off()
    expect_true(TRUE)
})

test_that("Test that plotDotPlot works with character class x variable", {
    temp <- tempfile("plot", fileext = ".png")
    png(filename = temp, width = 2750, height = 1600, res = 300)
    plotDotPlot(
        df = fedupPlot[1:10,],
        xVar = "set",
        yVar = "pathway",
        xLab = "-log10(qvalue)",
        fillVar = "sign",
        fillLab = "Genetic interaction",
        sizeVar = "fold_enrichment",
        sizeLab = "Fold enrichment"
    )
    dev.off()
    expect_true(TRUE)
})
