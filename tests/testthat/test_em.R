context("Enrichment map")

test_that("Test that writeFemap works", {
  data(testGene)
  data(backgroundGene)
  data(pathwaysGMT)
  fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
  results_file <- tempfile("fedup_res", fileext = ".txt")
  writeFemap(fedup_res, results_file)
  femap_res <- read.delim(results_file)

  expect_equal(nrow(fedup_res), nrow(femap_res))
  expect_true("status" %in% colnames(femap_res))
  expect_true(fedup_res[1,"status"] == "Enriched" && femap_res[1,"status"] == 1)
  expect_true(fedup_res[1436,"status"] == "Depleted" && femap_res[1436,"status"] == -1)
})
