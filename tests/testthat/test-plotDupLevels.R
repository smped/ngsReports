
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
  p <- plotDupLevels(fp)
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fp, TRUE)
  expect_true(is(p, "plotly"))
  p <- plotDupLevels(fpl, plotType = "bar")
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fpl, plotType = "heatmap")
  expect_true(is(p, "gg"))
  p <- plotDupLevels(fpl, TRUE, plotType = "bar")
  expect_true(is(p, "plotly"))
  p <- plotDupLevels(fpl, TRUE, plotType = "heatmap")
  expect_true(is(p, "plotly"))

})
