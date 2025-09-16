
test_that("Basic FastpData Plots work", {
  p <- plotInsertSize(fp)
  expect_true(is_ggplot(p))
  p <- plotInsertSize(fp, plotType = "cumulative")
  expect_true(is_ggplot(p))
  p <- plotInsertSize(fp, TRUE)
  expect_true(is(p, "plotly"))
  p <- plotInsertSize(fp, TRUE, plotType = "cumulative")
  expect_true(is(p, "plotly"))
})

test_that("Basic FastpDataList Plots work", {
  ## Heatmaps
  p <- plotInsertSize(fpl)
  expect_true(is_ggplot(p))
  expect_equal(p@labels$title, "Insert Size Distribution")
  expect_equal(
    lapply(p@mapping, rlang::as_label),
    list(
      x = "Insert Size", y = "Filename", fill = "Frequency", perc = "%",
      total = "Total"
    )
  )
  p <- plotInsertSize(fp, TRUE)
  expect_true(is(p, "plotly"))
  ## Lines
  p <- plotInsertSize(fpl, plotType = "l")
  expect_true(is_ggplot(p))
  expect_equal(
    lapply(p@mapping, rlang::as_label), list(x = "Insert Size", y = "Frequency")
  )
  expect_equal(
    lapply(p@layers$geom_line$mapping, rlang::as_label),
    list(colour = "Filename")
  )
  p <- plotInsertSize(fpl, plotType = "l", usePlotly = TRUE)
  expect_true(is(p, "plotly"))
  ## Cumulative
  p <- plotInsertSize(fpl, plotType = "c")
  expect_true(is_ggplot(p))
  expect_equal(
    lapply(c(p@mapping, p@layers$geom_line$mapping), rlang::as_label),
    list(x = "Insert Size", y = "Cumulative Frequency", colour = "Filename")
  )
  p <- plotInsertSize(fpl, plotType = "c", usePlotly = TRUE)
  expect_true(is(p, "plotly"))
})

test_that("Expected Errors appear", {
  expect_error(plotInsertSize(fp, plotTheme = ""))
  expect_error(plotInsertSize(fpl, plotType = "l", scaleColour = ""))
})
