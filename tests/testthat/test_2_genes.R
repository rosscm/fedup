context("(2) Genes")

test_that("Test that prepInput works", {
    genesFile <- system.file("extdata", "NIHMS1587165-supplement-1587165_Sup_Tab_2.txt", package = "fedup")
    genes <- read.delim(genesFile, h = TRUE, as.is = TRUE)

    b <- unique(genes[, "gene"])
    set1 <- unique(genes[which(genes$FASN_merge < -0.4), "gene"])
    set2 <- unique(genes[which(genes$FASN_merge > 0.4), "gene"])

    expect_error(prepInput(c("background", "negative"), b))
    expect_error(prepInput(NULL, b, set1))
    expect_error(prepInput(c("negative", "negative"), set1, set2))
    expect_error(prepInput(c("negative", "positive"), set1, set2))
    expect_true(is.list(prepInput(c("background", "negative", "positive"), b, set1, set2)))
})
