
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
