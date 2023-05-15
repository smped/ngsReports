## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)
fpl <- FastpDataList(system.file("extdata/fastp.json", package = "ngsReports"))

test_that("plotReadTotals works for FastQC objects", {
  p <- plotReadTotals(fdl)
  expect_true(is(p, "gg"))
  p <- plotReadTotals(fdl, duplicated = FALSE, scaleY = 1e3)
  expect_true(is(p, "gg"))
  p <- plotReadTotals(fdl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
})

test_that("plotReadtTotals works for Fastp objects", {
  p <- plotReadTotals(fpl)
  expect_true(is(p, "gg"))
  p <- plotReadTotals(fpl, status = FALSE, vertBars = TRUE)
  expect_true(is(p, "gg"))
  p <- plotReadTotals(fpl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
})
