packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastp.json", full.names = TRUE)
fp <- FastpData(fl)
fpl <- FastpDataList(path(fp))

test_that("Basic FastpData Plots work", {
  p <- plotInsertSize(fp)
  expect_true(is(p, "gg"))
  p <- plotInsertSize(fp, plotType = "cumulative")
  expect_true(is(p, "gg"))
  p <- plotInsertSize(fp, TRUE)
  expect_true(is(p, "plotly"))
  p <- plotInsertSize(fp, TRUE, plotType = "cumulative")
  expect_true(is(p, "plotly"))
})

test_that("Basic FastpDataList Plots work", {
  p <- plotInsertSize(fpl)
  expect_true(is(p, "gg"))
  p <- plotInsertSize(fp, TRUE)
  expect_true(is(p, "plotly"))
})

test_that("Expected Errors appear", {
  expect_error(plotInsertSize(fp, plotTheme = ""))
})
