#' @title Plot the Per Sequence GC Content
#'
#' @description Plot the Per Sequence GC Content for a set of FASTQC files
#'
#' @details
#' Makes plots for GC_Content.
#' When applied to a single FastqcFile or FastqcData object a simple line plot
#' will be drawn, with Theoretical GC content overlaid if desired
#'
#' For a FastqcFileList, or FastqcDataList either a line plot or heatmap can be
#' drawn
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList},
#' \code{FastqcData}, \code{FastqcDataList} or file path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using
#' ggplot. If \code{TRUE} plot will be rendered with plotly
#' @param counts \code{logical}. Plot the counts from each file if
#' \code{counts = TRUE}, otherwise frequencies will be plotted.
#' Ignored if calling the function on a FastqcDataList.
#' @param theoreticalGC \code{logical} default is \code{FALSE} to give the true
#' GC content, set to \code{TRUE} to normalize values of GC_Content by the
#' theoretical values using \code{\link{gcTheoretical}}. \code{species} must be
#' specified.
#' @param theoreticalType \code{character} Select type of data to normalize GC
#' content against. Accepts either "Genome" (default) or "Transcriptome".
#' @param GCobject an object of class GCTheoretical.
#'  Defaults to the gcTheoretical object supplied witht= the package
#' @param Fastafile a fasta file contains DNA sequences to generate theoretical
#' GC content
#' @param n number of simulated reads to generate theoretical GC content from
#' \code{Fastafile}
#' @param bp simulated read length to generate theoretical GC content from
#' \code{Fastafile}
#' @param species \code{character} if \code{gcTheory} is \code{TRUE} it must be
#' accompanied by a species. Species currently supported can be obtained using
#' \code{mData(gcTheoretical)}
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType Takes values "line" or "heatmap"
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for
#' pass, warning, and fail values in plot
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE}
#' then the dendrogram  will be displayed.
#' @param lineCols Colors for observed and theoretical GC lines in single plots
#' @param ... Used to pass various potting parameters to theme.
#'
#' @return A ggplot2 or plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot for a FastqcDataList
#' plotGcContent(fdl)
#'
#' # Plot a single FastqcData object
#' plotGcContent(fdl[[1]])
#'
#' # Plot GC content with theoretical GC content generated from a given
#' # fasta file
#' f <- "Athaliana.TAIR10.tRNA.fasta"
#' faFile <- system.file("extdata", f, package="ngsReports")
#' plotGcContent(fdl, Fastafile = faFile)
#'
#' @importFrom viridisLite inferno
#' @importFrom grDevices colorRampPalette
#' @importFrom stats hclust dist
#' @import ggplot2
#' @import fastqcTheoreticalGC
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocGenerics width
#' @importFrom utils read.table
#' @name plotGcContent
#' @rdname plotGcContent-methods
#' @export
setGeneric(
    "plotGcContent",
    function(
        x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
        theoreticalType = "Genome", species = "Hsapiens", GCobject, Fastafile,
        n = 1e+6, bp = 100, ...){
        standardGeneric("plotGcContent")
    }
)
#' @aliases plotGcContent,character
#' @rdname plotGcContent-methods
#' @export
setMethod(
    "plotGcContent",
    signature = "character",
    function(
        x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
        theoreticalType = "Genome", species = "Hsapiens", GCobject, Fastafile,
        n = 1e+6, bp = 100, ...
    ){
        x <- getFastqcData(x)
        plotGcContent(
            x, usePlotly, labels, theoreticalGC, theoreticalType, species,
            GCobject, Fastafile, n, bp, ...
        )
    }
)
#' @aliases plotGcContent,FastqcFile
#' @rdname plotGcContent-methods
#' @export
setMethod(
    "plotGcContent",
    signature = "FastqcFile",
    function(
        x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
        theoreticalType = "Genome", species = "Hsapiens", GCobject, Fastafile,
        n = 1e+6, bp = 100, ...
    ){
        x <- getFastqcData(x)
        plotGcContent(
            x, usePlotly, labels, theoreticalGC, theoreticalType, species,
            GCobject, Fastafile, n, bp, ...
        )
    }
)
#' @aliases plotGcContent,FastqcFileList
#' @rdname plotGcContent-methods
#' @export
setMethod(
    "plotGcContent",
    signature = "FastqcFileList",
    function(
        x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
        theoreticalType = "Genome", species = "Hsapiens", GCobject, Fastafile,
        n = 1e+6, bp = 100, ...
    ){
        x <- getFastqcData(x)
        plotGcContent(
            x, usePlotly, labels, theoreticalGC, theoreticalType,species,
            GCobject, Fastafile, n, bp, ...
        )
    }
)
#' @aliases plotGcContent,FastqcData
#' @rdname plotGcContent-methods
#' @export
setMethod(
    "plotGcContent",
    signature = "FastqcData",
    function(
        x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
        theoreticalType = "Genome", species = "Hsapiens", GCobject, Fastafile,
        n = 1e+6, bp = 100, counts = FALSE, lineCols = c("red", "blue"), ...
    ){

        df <- Per_sequence_GC_content(x)

        if (!length(df)) {
            gcPlot <- .emptyPlot("No Duplication Levels Module Detected")
            if (usePlotly) gcPlot <- ggplotly(gcPlot, tooltip = "")
            return(gcPlot)
        }

        df$Type <- "GC count per read"

        ## Get the correct y-axis label
        yLab <- c("Frequency", "Count")[counts + 1]

        ## Set the labels
        labels <- .makeLabels(df, labels, ...)
        df$Filename <- labels[df$Filename]

        ## Get any arguments for dotArgs that have been set manually
        dotArgs <- list(...)
        allowed <- names(formals(ggplot2::theme))
        keepArgs <- which(names(dotArgs) %in% allowed)
        userTheme <- c()
        if (length(keepArgs) > 0)
            userTheme <- do.call(theme, dotArgs[keepArgs])

        ## Tidy up the GC content variables
        if (missing(GCobject)) {
            GCobject <- ngsReports::gcTheoretical
        }

        if (theoreticalGC) {
            if (!missing(Fastafile)) {
                gcTheoryDF <- .gcFromFasta(Fastafile,n,bp)
                subTitle <-
                    paste("Theoretical Distribution based on file", Fastafile)
            }
            else{
                gcFun <- match.arg(
                    tolower(theoreticalType),
                    c("genomes","transcriptomes")
                )
                avail <- do.call(gcFun, list(object = GCobject))
                stopifnot(species %in% avail$Name)
                gcTheoryDF <-
                    getGC(GCobject, name = species, type = theoreticalType)
                names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
                subTitle <- paste(
                    "Theoretical Distribution based on the",
                    species,
                    theoreticalType
                )
            }
            gcTheoryDF$Type <- "Theoretical Distribution"
            gcTheoryDF$Filename <- "Theoretical Distribution"
            gcTheoryDF$Freq <- round(gcTheoryDF$Freq,4)
        }
        else{
            gcTheoryDF <- c()
            subTitle <- c()
        }

        xLab <- "GC Content (%)"

        if (!counts) {## If using frequencies (not counts)

            ## Summarise to frequencies & initialise the plot
            df$Freq <- df$Count / sum(df$Count)
            df <- df[c("Type", "GC_Content", "Freq")]
            df <- dplyr::bind_rows(df, gcTheoryDF)
            df$Type <- as.factor(df$Type)
            df$Freq <- round(df$Freq, 4)

            gcPlot <- ggplot(
                df,
                aes_string(x = "GC_Content", y = "Freq", colour = "Type")
            ) + geom_line()
        }
        else{

            df <- df[c("GC_Content","Type", "Count")]
            if (theoreticalGC) {
                gcTheoryDF$Count <- gcTheoryDF$Freq * sum(df$Count)
                gcTheoryDF <- gcTheoryDF[c("GC_Content","Type", "Count")]
                df <- dplyr::bind_rows(df, gcTheoryDF)
            }
            ## Initialise the plot using counts
            gcPlot <- ggplot(
                df,
                aes_string(x = "GC_Content",
                           y = "Count",
                           colour = "Type")
            ) + geom_line()
        }

        gcPlot <- gcPlot +
            scale_colour_manual(values = lineCols) +
            scale_x_continuous(
                breaks = seq(0, 100, by = 10), expand = c(0.02, 0)
            ) +
            labs(x = xLab, y = yLab, colour = c()) +
            ggtitle(label = labels,
                    subtitle = subTitle) +
            theme_bw() +
            theme(
                legend.position = c(1, 1),
                legend.justification = c(1, 1),
                legend.background = element_rect(
                    colour = "grey20",
                    size = 0.2
                ),
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)
            )
        if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

        if (usePlotly) {
            value <- c("Freq", "Count")[counts + 1]
            ttip <- c("GC_Content", value, "Type")
            gcPlot <- gcPlot + ggtitle(label = labels, subtitle = c())
            gcPlot <- suppressWarnings(
                suppressMessages(
                    ## Try using subplot with plotly_empty to align with
                    ## the heatmap in the app
                    plotly::ggplotly(gcPlot, tooltip = ttip)
                )
            )
            gcPlot <- plotly::layout(gcPlot, legend = list(x = 0.622, y = 1))
            gcPlot <- suppressMessages(
                plotly::subplot(
                    plotly::plotly_empty(),
                    gcPlot,
                    widths = c(0.14,0.86)
                )
            )
            gcPlot <- plotly::layout(
                gcPlot,
                xaxis2 = list(title = xLab),
                yaxis2 = list(title = yLab)
            )

        }

        ## Draw the plot
        gcPlot
    }
)
#' @aliases plotGcContent,FastqcDataList
#' @rdname plotGcContent-methods
#' @export
setMethod(
    "plotGcContent",
    signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
        theoreticalType = "Genome", species = "Hsapiens", GCobject, Fastafile,
        n=1e+6, bp=100, plotType = c("heatmap", "line"), pwfCols,
        cluster = FALSE, dendrogram = FALSE, ...
    ){

        df <- Per_sequence_GC_content(x)

        if (!length(df)) {
            gcPlot <- .emptyPlot("No Duplication Levels Module Detected")
            if (usePlotly) gcPlot <- ggplotly(gcPlot, tooltip = "")
            return(gcPlot)
        }

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(df, labels, ...)

        ## Always use frequencies not counts
        df <- lapply(split(df, f = df$Filename), function(x){
            x$Percent <- 100*x$Count / sum(x$Count)
            x
        })
        df <- dplyr::bind_rows(df)
        df <- df[c("Filename", "GC_Content", "Percent")]

        ## Get any arguments for dotArgs that have been set manually
        dotArgs <- list(...)
        if ("size" %in% names(dotArgs)) {
            lineWidth <- dotArgs$size
        }
        else{
            lineWidth <- c(line = 0.5, heatmap = 0.2)[plotType]
        }

        allowed <- names(formals(ggplot2::theme))
        keepArgs <- which(names(dotArgs) %in% allowed)
        userTheme <- c()
        if (length(keepArgs) > 0)
            userTheme <- do.call(theme, dotArgs[keepArgs])

        ## Initialise objects then fill if required
        gcTheoryDF <- c()
        subTitle <- c()
        if (theoreticalGC) {
            if (!missing(Fastafile)) {
                gcTheoryDF <- suppressMessages(
                    .gcFromFasta(Fastafile, n, bp)
                )
                subTitle <- paste(
                    "Theoretical Distribution based on",
                    basename(Fastafile))
            }
            else {
                ## Tidy up the GC content variables
                if (missing(GCobject)) {
                    GCobject <- gcTheoretical
                }
                gcFun <- match.arg(
                    tolower(theoreticalType),
                    c("genomes","transcriptomes")
                )
                avail <- do.call(gcFun, list(object = GCobject))
                species <- match.arg(species, avail$Name)
                gcTheoryDF <- getGC(
                    GCobject,
                    name = species,
                    type = theoreticalType
                )
                names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
                subTitle <- paste(
                    "Theoretical Distribution based on the",
                    species, theoreticalType
                )
            }
            gcTheoryDF$Filename <- "Theoretical Distribution"
            gcTheoryDF$Percent <- round(100*gcTheoryDF$Freq,4)
        }

        ## Check for valid plotType arguments
        plotType <- match.arg(plotType)

        xLab <- "GC Content (%)"

        if (plotType == "line") {
            df$Filename <- labels[df$Filename]
            ## Setup a palette with black as the first colour.
            ## Use the paired palette for easier visualisation of paired
            ## data
            n <- length(x)
            lineCols <- RColorBrewer::brewer.pal(min(12, n), "Paired")
            lineCols <- colorRampPalette(lineCols)(n)
            lineCols <- c("#000000", lineCols)

            df <- dplyr::bind_rows(gcTheoryDF, df)
            df$Filename <- factor(df$Filename, levels = unique(df$Filename))
            df$Percent <- round(df$Percent, 2)

            gcPlot <- ggplot(
                df,
                aes_string(
                    x = "GC_Content",
                    y = "Percent",
                    colour = "Filename"
                )
            ) +
                geom_line(size = lineWidth) +
                scale_x_continuous(
                    breaks = seq(0, 100, by = 10),
                    expand = c(0.02, 0)
                ) +
                scale_colour_manual(values = lineCols) +
                labs(x = xLab, y = "Reads (%)", colour = c()) +
                ggtitle(label = c(), subtitle = subTitle) +
                theme_bw() +
                theme(
                    plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5)
                )
            if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

            if (usePlotly) {

                gcPlot <- gcPlot +
                    labs(colour = "Filename") +
                    theme(legend.position = "none") +
                    ggtitle(NULL)
                gcPlot <- suppressWarnings(
                    suppressMessages(
                        plotly::ggplotly(
                            gcPlot,
                            tooltip = c("x", "y", "colour")
                        )
                    )
                )
            }
        }

        if (plotType == "heatmap") {

            if (theoreticalGC) {
                ## If using theoretical GC, just show the difference
                df <- lapply(split(df, df$Filename), function(x){
                    x$Percent <- x$Percent - unlist(gcTheoryDF$Percent)
                    x$Percent <- round(x$Percent, 2)
                    x
                })
                df <- dplyr::bind_rows(df)
                fillLab <- "Difference from\nTheoretical GC"
            }
            else{
                fillLab <- "Frequency"
            }

            if (dendrogram && !cluster) {
                message("cluster will be set to TRUE when dendrogram = TRUE")
                cluster <- TRUE
            }

            key <- names(labels)
            if (cluster) {
                clusterDend <-
                    .makeDendrogram(df, "Filename", "GC_Content", "Percent")
                key <- labels(clusterDend)
            }
            ## Now set everything as factors
            df$Filename <- factor(labels[df$Filename], levels = labels[key])


            ## Draw the heatmap
            gcPlot <-
                ggplot(
                    df,
                    aes_string("GC_Content","Filename", fill = "Percent")
                ) +
                geom_tile() +
                scale_x_continuous(expand = c(0, 0)) +
                theme(
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank()
                ) +
                labs(x = xLab, fill = fillLab) +
                scale_fill_gradient2(
                    low = inferno(1, begin = 0.4),
                    high = inferno(1, begin = 0.9),
                    midpoint = 0,
                    mid = inferno(1, begin = 0)
                )

            if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

            if (usePlotly) {

                gcPlot <- gcPlot +
                    theme(
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.position = "none"
                    )

                # Prepare the sideBar
                if (missing(pwfCols)) pwfCols <- ngsReports::pwf
                status <- getSummary(x)
                status <- subset(status, Category == "Per sequence GC content")
                status$Filename <- factor(
                    labels[status$Filename],
                    levels = levels(df$Filename)
                )
                status <- dplyr::right_join(
                    status,
                    unique(df["Filename"]),
                    by = "Filename"
                )
                sideBar <- .makeSidebar(status, key, pwfCols)

                ##plot dendrogram
                if (dendrogram) {
                    dx <- ggdendro::dendro_data(clusterDend)
                    dendro <- .renderDendro(dx$segments)
                }
                else{
                    dendro <- plotly::plotly_empty()
                }

                gcPlot <- suppressMessages(
                    suppressWarnings(
                        plotly::subplot(
                            dendro,
                            sideBar,
                            gcPlot,
                            widths = c(0.1,0.08,0.82),
                            margin = 0.001,
                            shareY = TRUE)
                    )
                )
                gcPlot <-
                    plotly::layout(gcPlot, xaxis3 = list(title = xLab))
            }
        }
        gcPlot
    }
)