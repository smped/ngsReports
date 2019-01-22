context("Check structure of the sideBar")

# This will give the correct structure as a stand alone process
packageDir <- system.file("extdata", package = "ngsReports")
fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)[1:2]
fdl <- getFastqcData(fileList)
status <- getSummary(fdl)
status <- subset(status, Category == "Basic Statistics")
key <- status$Filename
status$Filename <- gsub(".fastq", "", status$Filename)
sb <- .makeSidebar(status, key, pwf)

test_that("A plot can be drawn",{
    expect_equal(class(sb), c("plotly", "htmlwidget"))
})

test_that("Plot is as expected",{
    # This hash was found as the first two digits of the md5sum from
    # digest::digest(sb$x$data, "md5")
    expect_known_hash(sb$x$data, 91)
})

closeAllConnections()
