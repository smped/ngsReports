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

test_that("plotBaseQuals outputs correct objects", {

  p <- plotBaseQuals(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotBaseQuals(fdl)
  expect_true(is(p, "gg"))

  p <- plotBaseQuals(fdl, plotType = "boxplot")
  expect_true(is(p, "gg"))

  p <- plotBaseQuals(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotBaseQuals(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})


test_that("plotDupLevels outputs correct objects", {

  p <- plotDupLevels(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotDupLevels(fdl)
  expect_true(is(p, "gg"))

  p <- plotDupLevels(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotDupLevels(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotDupLevels(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotGcContent outputs correct objects", {

  p <- plotGcContent(fdl[[1]])
  expect_true(is(p, "gg"))

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

})

test_that("plotKmers outputs correct objects", {

  p <- plotKmers(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotKmers(fdl)
  expect_true(is(p, "gg"))

  p <- plotKmers(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotKmers(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotNContent outputs correct objects", {

  p <- plotNContent(fdl)
  expect_true(is(p, "gg"))

  p <- plotNContent(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotOverrep outputs correct objects", {

  p <- plotOverrep(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotOverrep(fdl)
  expect_true(is(p, "gg"))

  p <- plotOverrep(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotOverrep(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotReadTotals works", {
  p <- plotReadTotals(fdl)
  expect_true(is(p, "gg"))
  p <- plotReadTotals(fdl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
})

test_that("plotSeqContent outputs correct objects", {

  p <- plotSeqContent(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl)
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl, plotType = "residuals")
  expect_true(is(p, "gg"))

  p <- plotSeqContent(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotSeqContent(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})


test_that("plotSeqLengthDistn outputs correct objects", {

  p <- plotSeqLengthDistn(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl[[1]], plotType = "cdf")
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl)
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl, plotType = "cdf")
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotSeqLengthDistn(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotSeqQuals outputs correct objects", {

  p <- plotSeqQuals(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl)
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotSeqQuals(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotSummary works", {
  p <- plotSummary(fdl)
  expect_true(is(p, "gg"))
  p <- plotSummary(fdl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
})
