context("Test getSummary() for correct loading")

test_that("Example file is loaded correctly", {
  # Define the example file with known structure
  x <- system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports")
  df <- getSummary(x)
  expect_equal(names(df), c("Status", "Category", "Filename"))
  expect_equal(nrow(df), 12)
})


