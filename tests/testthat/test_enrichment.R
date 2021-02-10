context("Enrichment analysis")

test_that("Test that FEDUP analysis works", {
  data(testGene)
  data(backgroundGene)
  data(pathwaysGMT)
  fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)

  expect_equal(fedup_res[1, real_frac], 100.00000)
  expect_equal(fedup_res[1, qvalue], 1.567426e-186)
  expect_true("NKX2-5" %in% fedup_res[,real_gene][[1]])
  expect_true(!"OR11A1" %in% fedup_res[,real_gene][[1]])
  expect_equal(fedup_res[1437, real_frac], 0.0000000)
  expect_equal(fedup_res[1437, qvalue], 1.000000e+00)
  expect_false("NKX2-5" %in% fedup_res[,real_gene][[1437]])
  expect_true(!"OR11A1" %in% fedup_res[,real_gene][[1437]])
})
