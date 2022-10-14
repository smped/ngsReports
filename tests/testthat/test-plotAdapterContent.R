packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)

test_that("plotAdapterContent outputs correct object classes", {
  p <- plotAdapterContent(fdl, dendrogram = TRUE, cluster= TRUE)
  expect_true(is(p, "patchwork"))
  p <- plotAdapterContent(fdl, plotType = "line")
  expect_true(is(p, "gg"))
  p <- plotAdapterContent(fdl, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
  p <- plotAdapterContent(fdl, plotType = "line", usePlotly = TRUE)
  expect_true(is(p, "plotly"))
})

test_that("plotAdapterContent errors", {
  expect_error(plotAdapterContent(fdl, adapterType = ""), "'arg' should be")
})

test_that(".prepHeatmap errors correctly", {
  expect_error(.prepHeatmap(""), 'is\\(x, "gg"\\) is not TRUE')
  p <- ggplot()
  expect_error(
    .prepHeatmap(p, tibble(x = numeric())),
    'all\\(c\\("Filename", "Status"\\) %in% colnames\\(status\\)\\) is not TRUE'
  )
  expect_error(
    .prepHeatmap(
      p, tibble(Filename = character(), Status = character()), segments = ""
    ),
    'all\\(c\\("x", "y", "xend", "yend"\\) %in% colnames\\(segments\\)\\) is not TRUE'
  )
})
