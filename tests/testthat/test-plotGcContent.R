
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

  fl <- system.file("extdata", "fastp.json.gz", package = "ngsReports")
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
