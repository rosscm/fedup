context("Enrichment analysis")

test_that("Test that FEDUP analysis works", {
  data(testGene)
  data(backgroundGene)
  data(pathwaysGMT)
  fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)

  expect_equal(fedup_res[1, real_pathway_frac], 100.00000)
  expect_equal(fedup_res[1, expected_pathway_frac], 2.48042593)
  expect_equal(fedup_res[1, fdr], 0.000000e+00)
  expect_equal(fedup_res[8080, real_pathway_frac], 0)
  expect_equal(fedup_res[8080, expected_pathway_frac], 0.18791106)
  expect_equal(fedup_res[8080, fdr], 1.000000e+00)

  expect_true("OR11A1" %in% fedup_res[,real_pathway_gene][[1]])
  expect_false(!"RTP4" %in% fedup_res[,real_pathway_gene][[1]])
})
