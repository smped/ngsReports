context("Check Basic Structure of FastqcFile")
# TODO
test_that("FastqcFile Can Be Formed", {
  fl <- system.file("extdata", "ATTG_R1_fastqc.zip", package = "ngsReports")
  fqcFile <- FastqcFile(fl)
  expect_equal(length(fqcFile), 1)
})

closeAllConnections()
