test_that("Example file is loaded correctly", {
  # Define the example file with known structure
  df <- getSummary(fl[[1]])
  expect_equal(names(df), c("Status", "Category", "Filename"))
  expect_equal(nrow(df), 12)
})

closeAllConnections()
