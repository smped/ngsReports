test_that("FastpData objects behave as expected", {
  df <- summariseOverrep(fpl, min_count = 1e3)
  expect_true(nrow(df) == 1)
  expect_equal(
    colnames(df),
    c(
      "reads", "sequence", "n_samples", "Before_count_mean", "Before_count_sum",
      "Before_count_max", "After_count_mean", "After_count_sum", "After_count_max",
      "Before_rate_mean", "Before_rate_sum", "Before_rate_max", "After_rate_mean",
      "After_rate_sum", "After_rate_max"
    )
  )
})

test_that("FastQC objects behave as expected", {
  df <- summariseOverrep(fdl, min_count = 100, pattern = "No")
  expect_true(nrow(df) == 73)
  expect_equal(
    colnames(df),
    c(
      "Sequence", "Possible_Source", "n_samples", "Count_mean", "Count_sum",
      "Count_max", "Percentage_mean", "Percentage_sum", "Percentage_max"
    )
  )

  df <- dplyr::bind_rows(lapply(fdl, summariseOverrep))
  expect_equal(nrow(df), length(fdl))
  expect_equal(
    colnames(df),
    c(
      "Filename", "n_sequences", "Count_mean", "Count_sum", "Count_max",
      "Percentage_mean", "Percentage_sum", "Percentage_max"
    )
  )

  df <- summariseOverrep(fdl[[1]], min_count = Inf)
  expect_true(!nrow(df))

})
