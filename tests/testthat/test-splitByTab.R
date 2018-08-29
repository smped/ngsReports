context("Check correct behaviour for splitting lines into a data frame")

test_that("splitByTab errors correctly with uneven tab delimiters",{
  x <- c("A\tB", "C\tD\tE")
  expect_error(splitByTab(x))
})

test_that("splitByTab errors correctly with no tab delimiters",{
  x <- c("A", "B")
  expect_error(splitByTab(x))
})

test_that("splitByTab forms correct default columns",{
  x <- c("Col1\tCol2", "Val1\tVal2")
  df <- splitByTab(x)
  expect_true(
    setequal(names(df), c("Col1", "Col2")) &&
      df[["Col1"]] == "Val1" &&
      df[["Col2"]] == "Val2"
  )
})

test_that("splitByTab forms correct columns without names",{
  x <- c("Col1\tCol2", "Val1\tVal2")
  df <- splitByTab(x, firstRowToNames = FALSE)
  expect_true(
    setequal(names(df), c("V1", "V2")) &&
      setequal(df[["V1"]], c("Col1", "Val1")) &&
      setequal(df[["V2"]], c("Col2", "Val2"))
  )
})
