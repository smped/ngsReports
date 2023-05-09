## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)

test_that("plotDupLevels outputs correct objects from FastQC", {

  p <- plotDupLevels(fdl[[1]])
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotDupLevels(fdl)
  expect_true(is(p, "gg"))

  p <- plotDupLevels(fdl, plotType = "line")
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fdl, TRUE, plotType = "line")
  expect_true(is(p, "plotly"))

  p <- plotDupLevels(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotDupLevels(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotDupLevels fastp outputs are correct",{
  fl <- system.file("extdata", "fastp.json", package = "ngsReports")
  fp <- FastpData(fl)
  p <- plotDupLevels(fp)
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fp, TRUE)
  expect_true(is(p, "plotly"))
  fpl <- FastpDataList(fl)
  p <- plotDupLevels(fpl, plotType = "bar")
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fpl, plotType = "heatmap")
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fpl, TRUE, plotType = "bar")
  expect_true(is(p, "plotly"))
  p <- plotDupLevels(fpl, TRUE, plotType = "heatmap")
  expect_true(is(p, "plotly"))

})
