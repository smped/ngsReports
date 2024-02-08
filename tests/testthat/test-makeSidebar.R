status <- getSummary(fdl)
status <- subset(status, Category == "Basic Statistics")
key <- status$Filename
status$Filename <- gsub(".fastq", "", status$Filename)
sb <- .makeSidebar(status, key, pwf)

test_that("A plot can be drawn",{
    expect_true(is(sb, "plotly"))
    expect_true(is(sb, "htmlwidget"))
})

test_that(
    "Plot is as expected",
    {
        expected <- c(
            x1 = "0.5", x2 = "0.5", x3 = "1.5", x4 = "1.5", x5 = "0.5",
            y1 = "0.5", y2 = "1.5", y3 = "1.5", y4 = "0.5", y5 = "0.5",
            text = "Filename: ATTG_R1<br />Status: PASS",
            key1 = "ATTG_R1.fastq", key2 = "ATTG_R1.fastq", key3 = "ATTG_R1.fastq",
            key4 = "ATTG_R1.fastq", key5 = "ATTG_R1.fastq", type = "scatter",
            mode = "lines", line.width = "0.377952755905512", line.color = "transparent",
            line.dash = "solid", fill = "toself", fillcolor = "rgba(0,204,0,1)",
            hoveron = "fills", name = "PASS", legendgroup = "PASS", showlegend = "TRUE",
            xaxis = "x", yaxis = "y", hoverinfo = "text"
        )
        x <- unlist(sb$x$data[[1]])
        expect_equal(
            x[1:10],
            c(x1 = "0.5", x2 = "0.5", x3 = "1.5", x4 = "1.5", x5 = "0.5",
              y1 = "0.5", y2 = "1.5", y3 = "1.5", y4 = "0.5", y5 = "0.5")
        )
        expect_true(
            all(
                grepl("Status: PASS", x[["text"]]),
                grepl("Filename: ATTG_R1", x[["text"]])
            )
        )
        expect_equal(unique(x[grepl("key", names(x))]), "ATTG_R1.fastq")
        expect_equal(
            x[c(
                "type", "mode", "line.width", "line.color", "line.dash",
                "fill", "fillcolor", "hoveron", "name", "legendgroup",
                "showlegend", "xaxis", "yaxis", "hoverinfo"
                )],
            c(
                type = "scatter", mode = "lines",
                line.width = "0.377952755905512", line.color = "transparent",
                line.dash = "solid", fill = "toself",
                fillcolor = "rgba(0,204,0,1)", hoveron = "fills", name = "PASS",
                legendgroup = "PASS", showlegend = "TRUE",
                xaxis = "x", yaxis = "y", hoverinfo = "text"
            )
        )
    }
)

closeAllConnections()
