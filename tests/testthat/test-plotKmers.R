
test_that("plotKmers outputs correct objects", {

  p <- plotKmers(fdl[[1]])
  expect_true(is(p, "gg"))
  p <- plotKmers(fdl[[1]], TRUE)
  expect_true(is(p, "plotly"))
  p <- plotKmers(fdl)
  expect_true(is(p, "gg"))
  p <- plotKmers(fdl, TRUE)
  expect_true(is(p, "plotly"))

  p <- plotKmers(fdl, dendrogam = TRUE)
  expect_true(is(p, "patchwork"))
  p <- plotKmers(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("FastpData objects plot correctly", {
  p <- plotKmers(fp)
  expect_true(is(p, "gg"))
  p <- plotKmers(fp, readsBy = "mean")
  expect_true(is(p, "gg"))
  p <- plotKmers(fp, TRUE)
  expect_true(is(p, "plotly"))
})


test_that("Correct Errors are produced", {
  expect_message(plotKmers(fpl), "Method not.*")
  expect_error(plotKmers(fp, plotTheme = ""))
  expect_error(plotKmers(fp, scaleFill = scale_fill_grey())) # Discrete Scale
})
