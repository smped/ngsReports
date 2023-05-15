## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)
fp <- FastpData(system.file("extdata/fastp.json", package = "ngsReports"))


test_that("plotNContent outputs correct objects", {

  p <- plotNContent(fdl[[1]])
  expect_true(is(p, "gg"))
  p <- plotNContent(fdl)
  expect_true(is(p, "gg"))

  p <- plotNContent(fdl[[1]], dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))
  p <- plotNContent(fdl, dendrogam = TRUE, usePlotly = TRUE)
  expect_true(is(p, "plotly"))

})

test_that("plotNContent works for all variations of FastpData params", {
  p <- plotNContent(fp)
  expect_true(is(p, "gg"))
  p <- plotNContent(fp, TRUE)
  expect_true(is(p, "plotly"))
})

test_that("plotNContent works for all variations of FastpDataList params", {
  fpl <- FastpDataList(path(fp))
  p <- plotNContent(fpl)
  expect_true(is(p, "gg"))
  expect_message(
    plotNContent(fpl, dendrogram = TRUE), "Cannot cluster.+"
  )
  expect_message(
    plotNContent(fpl[1], reads = "read1", dendrogram = TRUE), "Cannot cluster one file.+"
  )
})
