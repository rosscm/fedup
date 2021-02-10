context("Enrichment dotplot")

test_that("Test that plotDotPlot works", {
  data(testGene)
  data(backgroundGene)
  data(pathwaysGMT)
  fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
  fedup_enr <- head(fedup_res[with(fedup_res, which(status == "Enriched")),], 10)
  fedup_dep <- head(fedup_res[with(fedup_res, which(status == "Depleted")),], 10)
  fedup_plot <- rbind(fedup_enr, fedup_dep)
  fedup_plot$log10fdr <- -log10(fedup_plot$fdr + 1e-10) # log10-transform FDR for plotting
  fedup_plot$pathway <- gsub("\\%.*", "", fedup_plot$pathway) # clean pathway names
  temp <- tempfile("plot", fileext = ".png")
  png(filename = temp, width = 2750, height = 1600, res = 300)
  plotDotPlot(df = fedup_plot,
              x_var = "log10fdr",
              y_var = "pathway",
              x_lab = "-log10(FDR)",
              fill_var = "status",
              fill_lab = "Enrichment status",
              size_var = "fold_enrichment",
              size_lab = "Fold enrichment")
   dev.off()
   expect_true(TRUE)
})
