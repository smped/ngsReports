## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)

test_that("plotOverrep outputs correct objects", {

  p <- plotOverrep(fdl[[1]])
  expect_true(is(p, "gg"))
  p <- plotOverrep(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotOverrep(fdl, showPwf = FALSE)
  expect_true(is(p, "gg"))

  p <- plotOverrep(fdl)
  expect_true(is(p, "patchwork"))

  p <- plotOverrep(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotOverrep(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})



test_that("plotSeqLengthDistn outputs correct objects", {

  p <- plotSeqLengthDistn(path(fdl[[1]]))
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl[[1]], counts = FALSE)
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl[[1]], plotType = "cdf")
  expect_true(is(p, "gg"))

  p <- plotSeqLengthDistn(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotSeqLengthDistn(fdl[[1]], TRUE, plotType = "cdf")
  expect_true(is(p, "plotly"))

  p <- plotSeqLengthDistn(fdl[[1]], TRUE, plotType = "line")
  expect_true(is(p, "plotly"))

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

  p <- plotSeqQuals(path(fdl[[1]]))
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl[[1]], counts = TRUE)
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotSeqQuals(fdl)
  expect_true(is(p, "patchwork"))

  p <- plotSeqQuals(fdl, showPwf = FALSE)
  expect_true(is(p, "gg"))

  p <- plotSeqQuals(fdl, plotType = "line")
  expect_true(is(p, "gg"))
  p <- plotSeqQuals(fdl, TRUE, plotType = "line")
  expect_true(is(p, "plotly"))

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


test_that(".updateThemeFromDots works", {
  p <- ggplot()
  expect_equal(p$theme, list())
  expect_equal(.updateThemeFromDots(p)$theme, list())
  p <- .updateThemeFromDots(p, a = "a")
  expect_equal(.updateThemeFromDots(p)$theme, list())

  p <- .updateThemeFromDots(p, plot.title = element_blank())
  expect_true(is(p$theme$plot.title, "element_blank"))
})
