context("Test the isCompressed() behaves correctly")

test_that("TRUE statement from a zipped file",{
  fl <- system.file("extdata", "ATTG_R1_fastqc.zip", package = "ngsReports")
  expect_true(isCompressed(fl))
})

test_that("FALSE statement from a directory",{
  d <- system.file("extdata", package = "ngsReports")
  expect_false(isCompressed(d))
})

test_that("FALSE statement from a gz file using zip",{
  fl <- system.file("extdata", "errorTestingFiles", "hello.txt.gz", package = "ngsReports")
  expect_false(isCompressed(fl, type = "zip"))
})

test_that("TRUE statement from a gz file using zip",{
  fl <- system.file("extdata", "errorTestingFiles", "hello.txt.gz", package = "ngsReports")
  expect_true(isCompressed(fl, type = "gz"))
})
