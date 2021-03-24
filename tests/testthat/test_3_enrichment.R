context("(3) Analysis")

data(pathwaysGMT)
test_that("Test that fedup only works with proper inputs", {
    data(geneSingle)
    expect_true(is.list(pathwaysGMT))
    expect_error(runFedup(NULL, pathwaysGMT))
    expect_error(runFedup(geneSingle, NULL))

    geneSingle_noBgrd <- geneSingle
    geneSingle_noBgrd[["background"]] <- NULL
    expect_error(runFedup(geneSingle_noBgrd, pathwaysGMT))

    geneSingle_onlyBgrd <- geneSingle
    geneSingle_onlyBgrd[["FASN_negative"]] <- NULL
    expect_error(runFedup(geneSingle_onlyBgrd, pathwaysGMT))

    geneSingle_twoBgrd <- geneSingle
    geneSingle_twoBgrd[["background2"]] <- geneSingle[["background"]]
    names(geneSingle_twoBgrd)[3] <- "background"
    expect_error(runFedup(geneSingle_twoBgrd, pathwaysGMT))

    geneSingle_switch <- geneSingle
    names(geneSingle_switch)[1] <- "negative"
    names(geneSingle_switch)[2] <- "background"
    expect_error(runFedup(geneSingle_switch, pathwaysGMT))

    geneSingle_overlap <- geneSingle
    geneSingle_overlap[["FASN_positive"]] <- geneSingle[["FASN_negative"]]
    expect_warning(runFedup(geneSingle_overlap, pathwaysGMT))
})

test_that("Test that fedup analysis works with single test set", {
    data(geneSingle)
    fedupRes <- runFedup(geneSingle, pathwaysGMT)
    expect_equal(length(which(names(geneSingle) != "background")), length(fedupRes))
    expect_equal(fedupRes[[1]][1, status], "enriched")
    expect_equal(fedupRes[[1]][1, qvalue], 2.294321e-09)
    expect_true("ALG3" %in% fedupRes[[1]][, real_gene][[1]])
    expect_equal(fedupRes[[1]][1437, status], "depleted")
    expect_equal(fedupRes[[1]][1437, qvalue], 1.000000e+00)
    expect_false("ALG3" %in% fedupRes[[1]][, real_gene][[1437]])
})

test_that("Test that fedup analysis works with two test sets", {
    data(geneDouble)
    fedupRes <- runFedup(geneDouble, pathwaysGMT)
    expect_equal(length(which(names(geneDouble) != "background")), length(fedupRes))
    expect_equal(fedupRes[[1]][1, status], "enriched")
    expect_equal(fedupRes[[1]][1, qvalue], 2.294321e-09)
    expect_true("MOGS" %in% fedupRes[[1]][, real_gene][[1]])
    expect_equal(fedupRes[[1]][1437, status], "depleted")
    expect_equal(fedupRes[[1]][1437, qvalue], 1.000000e+00)
    expect_false("MOGS" %in% fedupRes[[1]][, real_gene][[1437]])
    expect_equal(fedupRes[[2]][1, status], "enriched")
    expect_equal(fedupRes[[2]][1, qvalue], 8.562503e-15)
    expect_true("EIF3A" %in% fedupRes[[2]][, real_gene][[1]])
    expect_equal(fedupRes[[2]][1437, status], "depleted")
    expect_equal(fedupRes[[2]][1437, qvalue], 1.000000e+00)
    expect_false("EIF3A" %in% fedupRes[[2]][, real_gene][[1437]])
})

test_that("Test that fedup analysis works with multiple test sets", {
    data(geneMulti)
    expect_warning(fedupRes <- runFedup(geneMulti, pathwaysGMT))
    expect_equal(length(which(names(geneMulti) != "background")), length(fedupRes))
    expect_true("MOGS" %in% fedupRes[[1]][, real_gene][[1]])
    expect_true("EIF3A" %in% fedupRes[[2]][, real_gene][[1]])
    expect_true("MOGS" %in% fedupRes[[3]][, real_gene][[1]])
    expect_true("RAD51AP1" %in% fedupRes[[4]][, real_gene][[1]])
    expect_true("ACAT2" %in% fedupRes[[5]][, real_gene][[1]])
    expect_true("SMARCD1" %in% fedupRes[[6]][, real_gene][[1]])
    expect_true("SUZ12" %in% fedupRes[[7]][, real_gene][[1]])
    expect_true("CARS2" %in% fedupRes[[8]][, real_gene][[1]])
    expect_true("MBTPS2" %in% fedupRes[[9]][, real_gene][[1]])
    expect_true("PEX19" %in% fedupRes[[10]][, real_gene][[1]])
    expect_true("TIMM23" %in% fedupRes[[11]][, real_gene][[1]])
    expect_true("RPN1" %in% fedupRes[[12]][, real_gene][[1]])
})
