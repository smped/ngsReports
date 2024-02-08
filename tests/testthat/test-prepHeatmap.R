
test_that(".prepHeatmap errors correctly", {
  expect_error(.prepHeatmap(""), 'is\\(x, "gg"\\) is not TRUE')
  p <- ggplot()
  # expect_error(
  #   .prepHeatmap(p, tibble(x = numeric())),
  #   'all\\(c\\("Filename", "Status"\\) %in% colnames\\(status\\)\\) is not TRUE'
  # )
  expect_error(
    .prepHeatmap(
      p, tibble(Filename = character(), Status = character()), segments = ""
    ),
    'all\\(c\\("x", "y", "xend", "yend"\\) %in% colnames\\(segments\\)\\) is not TRUE'
  )
})
