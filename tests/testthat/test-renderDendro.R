context("Check structure of the dendrogram")

# This will give the correct structure as a stand alone process
packageDir <- system.file("extdata", package = "ngsReports")
fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1:2]
fdl <- getFastqcData(fileList)

#test dataset
df <- getModule(fdl, "Sequence_Length_Distribution")
cols <- c("Filename", "Lower", "Count")

#cluster the data
clusterDend <- .makeDendro(df[cols], "Filename","Lower", "Count")

#build and render the dendrogram
dx <- ggdendro::dendro_data(clusterDend)
dendro <- .renderDendro(dx$segments)

#test that the plot is plotly class
test_that("A plot can be drawn",{
  expect_equal(class(dendro), c("plotly", "htmlwidget"))
})

test_that("Plot is as expected",{
    expected <- c(x1 = "0", x2 = "0", x3 = NA, x4 = "0", x5 = "0", x6 = NA, x7 = "0",
                  x8 = "0", x9 = NA, x10 = "0", x11 = "0", y1 = "1.5", y2 = "1",
                  y3 = NA, y4 = "1", y5 = "1", y6 = NA, y7 = "1.5", y8 = "2", y9 = NA,
                  y10 = "2", y11 = "2", text = "", type = "scatter", mode = "lines",
                  line.width = "1.88976377952756", line.color = "rgba(0,0,0,1)",
                  line.dash = "solid", hoveron = "points", showlegend = "FALSE",
                  xaxis = "x", yaxis = "y", hoverinfo = "text", name = "")
  expect_equal(unlist(dendro$x$data[[1]]), expected)
})

closeAllConnections()
