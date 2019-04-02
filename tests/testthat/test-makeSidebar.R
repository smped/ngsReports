context("Check structure of the sideBar")

# This will give the correct structure as a stand alone process
packageDir <- system.file("extdata", package = "ngsReports")
fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1:2]
fdl <- FastqcDataList(fileList)
status <- getSummary(fdl)
status <- subset(status, Category == "Basic Statistics")
key <- status$Filename
status$Filename <- gsub(".fastq", "", status$Filename)
sb <- .makeSidebar(status, key, pwf)

test_that("A plot can be drawn",{
    expect_equal(class(sb), c("plotly", "htmlwidget"))
})

test_that("Plot is as expected",{
    expected <- c(x1 = "0.5", x2 = "0.5", x3 = "1.5", x4 = "1.5", x5 = "0.5",
                  y1 = "0.5", y2 = "1.5", y3 = "1.5", y4 = "0.5", y5 = "0.5", text = "Status: PASS<br />Filename: ATTG_R1",
                  key1 = "ATTG_R1.fastq", key2 = "ATTG_R1.fastq", key3 = "ATTG_R1.fastq",
                  key4 = "ATTG_R1.fastq", key5 = "ATTG_R1.fastq", type = "scatter",
                  mode = "lines", line.width = "0.377952755905512", line.color = "transparent",
                  line.dash = "solid", fill = "toself", fillcolor = "rgba(0,204,0,1)",
                  hoveron = "fills", name = "PASS", legendgroup = "PASS", showlegend = "TRUE",
                  xaxis = "x", yaxis = "y", hoverinfo = "text")
    expect_equal(unlist(sb$x$data[[1]]), expected)
})

closeAllConnections()
