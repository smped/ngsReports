## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)
fp <- FastpData(system.file("extdata/fastp.json.gz", package = "ngsReports"))

test_that("plotBaseQuals (FastQC) outputs correct objects", {

  p <- plotBaseQuals(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotBaseQuals(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotBaseQuals(fdl)
  expect_true(is(p, "gg"))

  p <- plotBaseQuals(fdl, plotType = "boxplot")
  expect_true(is(p, "gg"))

  p <- plotBaseQuals(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotBaseQuals(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

  p <- plotBaseQuals(fdl, TRUE, plotType = "boxplot")
  expect_true(is(p, "plotly"))

})

test_that("plotBaseQuals (Fastp) output correct objects",{
  p <- plotBaseQuals(fp)
  expect_true(is(p, "gg"))
  p <- plotBaseQuals(fp, TRUE)
  expect_true(is(p, "plotly"))

  fpl <- FastpDataList(path(fp))
  p <- plotBaseQuals(fpl)
  expect_true(is(p, "gg"))
  p <- plotBaseQuals(fpl, TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotBaseQuals errors correctly", {
  expect_error(plotBaseQuals(fp, plotTheme = ""))
  expect_error(plotBaseQuals(fp, scaleColour = ""))
})
