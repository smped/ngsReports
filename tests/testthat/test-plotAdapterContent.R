## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)

test_that("plotAdapterContent outputs correct object classes", {

  p <- plotAdapterContent(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotAdapterContent(fdl, dendrogram = TRUE, cluster= TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotAdapterContent(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotAdapterContent(fdl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

  p <- plotAdapterContent(fdl, plotType = "line", usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})


test_that("FastpData plots correctly", {

  fl <- system.file("extdata", "fastp.json", package = "ngsReports")
  fp <- FastpData(fl)
  p <- plotAdapterContent(fp)
  expect_true(is(p, "gg"))
  p <- plotAdapterContent(fp, TRUE)
  expect_true(is(p, "plotly"))

})
