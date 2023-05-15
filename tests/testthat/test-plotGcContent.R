## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)


test_that("plotGcContent outputs correct objects", {

  p <- plotGcContent(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotGcContent(fdl[[1]], usePlotly = TRUE)
  expect_true(is(p, "plotly"))

  p <- plotGcContent(fdl)
  expect_true(is(p, "gg"))

  p <- plotGcContent(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotGcContent(fdl, plotType = "cdf")
  expect_true(is(p, "gg"))

  p <- plotGcContent(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotGcContent(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

  p <- plotGcContent(fdl, plotType = "cdf", usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotGcContent works for Fastp structures", {

  fl <- system.file("extdata", "fastp.json", package = "ngsReports")
  fpl <- FastpDataList(fl)

  p <- plotGcContent(fpl[[1]])
  expect_true(is(p, "gg"))
  p <- plotGcContent(fpl)
  expect_true(is(p, "gg"))

  p <- plotGcContent(fpl[[1]], usePlotly = TRUE)
  expect_true(is(p, "plotly"))
  p <- plotGcContent(fpl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})
