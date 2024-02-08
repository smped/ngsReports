
test_that("plotAdapterContent outputs correct object classes", {

  p <- plotAdapterContent(fdl[[1]])
  expect_true(is(p, "gg"))

  p <- plotAdapterContent(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))

  p <- plotAdapterContent(fdl, dendrogram = TRUE, cluster = TRUE)
  expect_true(is(p, "patchwork"))

  p <- plotAdapterContent(fdl, plotType = "line")
  expect_true(is(p, "gg"))

  p <- plotAdapterContent(fdl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
  expect_true(length(p$x$data) == 8)
  p <- plotAdapterContent(fdl, usePlotly = TRUE, showPwf = FALSE)
  expect_true(length(p$x$data) == 1)

  p <- plotAdapterContent(fdl, plotType = "line", usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})


test_that("FastpData plots correctly", {

  fl <- system.file("extdata", "fastp.json.gz", package = "ngsReports")
  fp <- FastpData(fl)
  p <- plotAdapterContent(fp)
  expect_true(is(p, "gg"))
  p <- plotAdapterContent(fp, TRUE)
  expect_true(is(p, "plotly"))
  expect_true(length(p$x$data) == 2)

})

test_that("FastpDataList plots correctly", {

  fl <- system.file("extdata", "fastp.json.gz", package = "ngsReports")
  fp <- FastpData(fl)
  fpl <- FastpDataList(path(fp))
  p <- plotAdapterContent(fpl)
  expect_true(is(p, "gg"))
  p <- plotAdapterContent(fpl, TRUE, dendrogram = TRUE, showPwf = TRUE)
  expect_true(is(p, "plotly"))
  expect_true(length(p$x$data) == 5)

})

test_that("Errors appear correct", {

  expect_message(plotAdapterContent(NULL), "Method not implemented.+")
  p <- plotAdapterContent(fdl, adapterType = "Next")
  expect_true(is(p$data, "waiver"))

})
