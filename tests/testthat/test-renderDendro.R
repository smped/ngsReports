context("Check structure of the dendrogram")

# This will give the correct structure as a stand alone process
packageDir <- system.file("extdata", package = "ngsReports")
fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)[1:2]
fdl <- getFastqcData(fileList)

#test dataset

df <- Sequence_Length_Distribution(fdl)
cols <- c("Filename", "Lower", "Count")

#cluster the data

clusterDend <-
  .makeDendrogram(df[cols], 
                  "Filename","Lower", "Count")


#build and render the dendrogram
dx <- ggdendro::dendro_data(clusterDend)
dendro <- .renderDendro(dx$segments)


#test that the plot is plotly class
test_that("A plot can be drawn",{
  expect_equal(class(dendro), c("plotly", "htmlwidget"))
})

test_that("Plot is as expected",{
  # This has was found as the first two digits of the md5sum from
  # digest::digest(dendro$x$data, "md5")
  expect_known_hash(dendro$x$data, "c7")
})

closeAllConnections()
