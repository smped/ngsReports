test_that(".makeDendro errors correctly with missing columns",{
    df <- data.frame()
    expect_error(.makeDendro(df))
})

test_that(".makeDendro errors correctly with incorrect columns",{
    df <- data.frame()
    expect_error(.makeDendro(df, "X", "Y", "Z"))
})

df <- data.frame(
  Filename = rep(c("A", "B"), times = 2), Position = rep(c(1, 2), each = 2),
  Value = 0
)
clus <- .makeDendro(df, rowVal = "Filename", colVal = "Position", value = "Value")

test_that(".makeDendro produces a correct dendrogram",{
    expect_equal(class(clus), "dendrogram")
    expect_equal(labels(clus), c("A", "B"))
    expect_equal(nobs(clus), 2)
    expect_equal(unlist(dendrapply(clus, function(x){attributes(x)$height})), c(0, 0))
})
