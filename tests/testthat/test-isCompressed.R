test_that("TRUE statement from a zipped file",{
  expect_true(isCompressed(fl[[1]]))
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

closeAllConnections()
