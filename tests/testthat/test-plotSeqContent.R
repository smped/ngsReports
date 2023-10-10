## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)
fp <- FastpData(system.file("extdata/fastp.json.gz", package = "ngsReports"))
fpl <- FastpDataList(path(fp))

test_that("plotSeqContent outputs correct objects", {

  p <- plotSeqContent(fdl[[1]])
  expect_true(is(p, "gg"))
  p <- plotSeqContent(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotSeqContent(fdl)
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl, plotType = "residuals")
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotSeqContent(fdl[1:4], dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotSeqContent(FastpData) works", {
  p <- plotSeqContent(fp)
  expect_true(is(p, "gg"))
  p <- plotSeqContent(fp, TRUE)
  expect_true(is(p, "plotly"))
})

test_that("plotSeqContent(FastpDataList) works", {
  p <- plotSeqContent(fpl)
  expect_true(is(p, "gg"))
  p <- plotSeqContent(fpl, showPwf = TRUE)
  expect_true(is(p, "patchwork"))
  p <- plotSeqContent(fpl, readsBy = "line", plotType = "line")
  expect_true(is(p, "gg"))
  p <- plotSeqContent(fpl, readsBy = "line", plotType = "resid")
  expect_true(is(p, "gg"))
  p <- plotSeqContent(fpl, TRUE)
  expect_true(is(p, "plotly"))
})

test_that("Errors are caught", {

  expect_error(
    plotSeqContent(fp, readsBy = "line", moduleBy = "line"),
    "Cannot set the same plotting parameter to both reads and module"
  )
  expect_error(
    plotSeqContent(fpl, plotType = "line"),
    "Cannot set reads and module to the same parameter for line plots"
  )
  expect_message(
    plotSeqContent(fpl, dendrogram = TRUE, reads = "read1", bases = "A"),
    "Cannot cluster one file. Ignoring cluster and dendgrogram"
  )

})
