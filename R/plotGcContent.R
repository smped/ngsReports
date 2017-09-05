#' @title Plot the Per Sequence GC Content
#'
#' @description Plot the Per Sequence GC Content for a set of FASTQC files
#'
#' @details
#' Makes plots for GC_Content.
#' When applied to a single FastqcFile or FastqcData object a simple line plot will be drawn,
#' with Theoretical GC content overlaid if desired
#'
#' For a FastqcFileList, or FastqcDataList either a line plot or heatmap can be drawn
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param theoreticalGC \code{logical} default is \code{FALSE} to give the true GC content%, set to \code{TRUE} to normalize
#' values of GC_Content by the theoretical values using \code{\link{gcTheoretical}}. \code{species} must be specified.
#' @param theoreticalType \code{character} Select type of data to normalize GC content% agianst accepts either "Genome" or
#' "Transcriptome". Default is "Genome"
#' @param GCobject an object of class GCTheoretical.
#'  Defaults to the gcTheoretical object supplied witht= the package
#' @param species \code{character} if \code{gcTheory} is \code{TRUE} its must be accompanied by a species
#' Species currently supported can be obtained using \code{mData(gcTheoretical)}
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType Takes values "line" or "heatmap"
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param lineCols Colors for observed and theoretical GC lines in single plots
#' @param ... Used to pass various potting parameters to theme.
#'
#' @return A ggplot2 or plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot for a FastqcDataList
#' plotGcContent(fdl)
#'
#' # Plot a single FastqcData object
#' plotGcContent(fdl[[1]])
#'
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom viridisLite inferno
#' @importFrom grDevices colorRampPalette
#'
#' @name plotGcContent
#' @rdname plotGcContent-methods
#' @export
setGeneric("plotGcContent",function(x, usePlotly = FALSE, labels, theoreticalGC = TRUE, theoreticalType = "Genome",
                                    species = "Hsapiens", GCobject, ...){standardGeneric("plotGcContent")})
#' @aliases plotGcContent,character
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "character",
          function(x, usePlotly = FALSE, labels, theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, ...){
            x <- getFastqcData(x)
            plotGcContent(x, usePlotly, labels = labels, theoreticalGC = theoreticalGC,
                          theoreticalType = theoreticalType, species = species, GCobject = GCobject, ...)
          }
)
#' @aliases plotGcContent,FastqcFile
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, labels, theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, ...){
            x <- getFastqcData(x)
            plotGcContent(x, usePlotly, labels = labels, theoreticalGC = theoreticalGC,
                          theoreticalType = theoreticalType, species = species, GCobject = GCobject, ...)
          }
)
#' @aliases plotGcContent,FastqcFileList
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, labels, theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, ...){
            x <- getFastqcData(x)
            plotGcContent(x, usePlotly, labels = labels, theoreticalGC = theoreticalGC,
                          theoreticalType = theoreticalType, species = species, GCobject = GCobject, ...)

          }
)
#' @aliases plotGcContent,FastqcData
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, counts = FALSE, lineCols = c("red", "blue"),
                   ...){

            df <- tryCatch(Per_sequence_GC_content(x))
            df$Type <- "GC count per read"

            # Get the correct y-axis label
            ylab <- c("Frequency", "Count")[counts + 1]

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                                  names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
            labels <- labels[df$Filename]

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Tidy up the GC content variables
            if(missing(GCobject)){
              GCobject <- ngsReports::gcTheoretical
            }

            if (theoreticalGC){
              gcFun <- tryCatch(match.arg(tolower(theoreticalType), c("genomes","transcriptomes")))
              avail <- do.call(gcFun, list(object = GCobject))
              stopifnot(species %in% avail$Name)
              gcTheoryDF <- getGC(GCobject, name = species, type = theoreticalType)
              names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
              gcTheoryDF$Type <- "Theoretical Distribution"
              subTitle <- paste("Theoretical Distribution based on the", species, theoreticalType)
            }
            else{
              gcTheoryDF <- c()
              subTitle <- c()
            }

            if (!counts){
              # If using frequencies
              # Summarise to frequencies & initialise the plot
              df$Freq <- df$Count / sum(df$Count)
              df <- df[c("Type", "GC_Content", "Freq")]
              df <- dplyr::bind_rows(df, gcTheoryDF)
              df$Type <- as.factor(df$Type)
              df$Freq <- round(df$Freq, 4)

              gcPlot <- ggplot(df, aes_string(x = "GC_Content", y = "Freq", colour = "Type")) +
                geom_line()
            }
            else{

              df <- df[c("GC_Content","Type", "Count")]
              if (theoreticalGC){
                gcTheoryDF$Count <- gcTheoryDF$Freq * sum(df$Count)
                gcTheoryDF <- gcTheoryDF[c("GC_Content","Type", "Count")]
                df <- dplyr::bind_rows(df, gcTheoryDF)
              }
              # Initialise the plot using counts
              gcPlot <- ggplot(df, aes_string(x = "GC_Content", y = "Count", colour = "Type")) +
                geom_line()

            }

            gcPlot <- gcPlot +
              scale_colour_manual(values = lineCols) +
              scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0.02, 0)) +
              labs(x = "GC Content",
                   y = ylab,
                   colour = c()) +
              ggtitle(label = labels,
                      subtitle = subTitle) +
              theme_bw() +
              theme(legend.position = c(1, 1),
                    legend.justification = c(1, 1),
                    legend.background = element_rect(colour = "grey20", size = 0.2),
                    plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
            if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

            if(usePlotly){
              value <- c("Freq", "Count")[counts + 1]
              gcPlot <- gcPlot +
                theme(legend.position = "none") +
                ggtitle(label = labels, subtitle = c())
              gcPlot <- suppressWarnings(
                suppressMessages(
                  plotly::ggplotly(gcPlot, tooltip = c("GC_Content", value, "Type"))
                )
              )
            }

            # Draw the plot
            gcPlot
          }
)
#' @aliases plotGcContent,FastqcDataList
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject,  plotType = "heatmap", pwfCols,
                   clusterNames = FALSE, dendrogram = TRUE, ...){

            df <- tryCatch(Per_sequence_GC_content(x))
            df$Type <- "GC count per read"

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                                  names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")


            # Always use frequencies
            df <- lapply(split(df, f = df$Filename), function(x){
              x$Freq <- x$Count / sum(x$Count)
              x
            })
            df <- dplyr::bind_rows(df)
            df <- df[c("Filename", "GC_Content", "Freq", "Type")]

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            if ("size" %in% names(dotArgs)){
              lineWidth <- dotArgs$size
            }
            else{
              lineWidth <- c(line = 0.5, heatmap = 0.2)[plotType]
            }
            if ("colour" %in% names(dotArgs) || "color" %in% names(dotArgs)){
              i <- which(names(dotArgs) %in% c("colour", "color"))
              lineCol <- dotArgs[[i]]
            }
            else{
              lineCol <- "grey20"
            }
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Tidy up the GC content variables
            if(missing(GCobject)){
              GCobject <- ngsReports::gcTheoretical
            }

            if (theoreticalGC){
              gcFun <- tryCatch(match.arg(tolower(theoreticalType), c("genomes","transcriptomes")))
              avail <- do.call(gcFun, list(object = GCobject))
              stopifnot(species %in% avail$Name)
              gcTheoryDF <- getGC(GCobject, name = species, type = theoreticalType)
              names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
              gcTheoryDF$Filename <- "Theoretical Distribution"
              gcTheoryDF$Type <- "Theoretical Distribution"
              subTitle <- paste("Theoretical Distribution based on the", species, theoreticalType)
              gcTheoryDF$Freq <- round(gcTheoryDF$Freq, 4)
            }
            else{
              gcTheoryDF <- c()
              subTitle <- c()
            }

            # Check for valid plotType arguments
            plotType <- match.arg(plotType, c("line", "heatmap"))

            if (plotType == "line"){
              df$Filename <- labels[df$Filename]
              # Setup a palette with black as the first colour.
              # Use the paired palette for easier visualisation of paired data
              n <- length(x)
              lineCols <- colorRampPalette(RColorBrewer::brewer.pal(min(12, n), "Paired"))(n)
              lineCols <- c("#000000", lineCols)

              df <- dplyr::bind_rows(gcTheoryDF, df)
              df$Filename <- factor(df$Filename, levels = unique(df$Filename))
              df$Freq <- round(df$Freq, 4)

              gcPlot <- ggplot(df, aes_string(x = "GC_Content", y = "Freq", colour = "Filename")) +
                geom_line(size = lineWidth) +
                scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0.02, 0)) +
                scale_colour_manual(values = lineCols) +
                labs(x = "GC Content",
                     y = "Frequency",
                     colour = c()) +
                ggtitle(label = c(),
                        subtitle = subTitle) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
              if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

              if(usePlotly){

                gcPlot <- gcPlot +
                  theme(legend.position = "none") +
                  ggtitle(c())
                gcPlot <- suppressWarnings(
                  suppressMessages(
                    plotly::ggplotly(gcPlot, tooltip = c("x", "y", "color", "Filename"))
                  )
                )
              }
            }

            if (plotType == "heatmap"){

              if (theoreticalGC){
                df <- lapply(split(df, df$Filename), function(x){
                  gcTheoryDF <- getGC(GCobject, name = species, type = theoreticalType)
                  x$Freq <- x$Freq - unlist(gcTheoryDF[species])
                  x
                }) %>%
                  dplyr::bind_rows()
                fillLab <- "Difference from\nTheoretical GC"
              }
              else{
                fillLab <- "Frequency"
              }

              if(clusterNames){
                # Grab the main columns & cast from long to wide
                mat <- reshape2::acast(df[c("Filename", "GC_Content", "Freq")], Filename ~ GC_Content, value.var = "Freq")
                mat[is.na(mat)] <- 0
                clus <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
                row.ord <- order.dendrogram(clus)
                mat <- mat[row.ord,]
                # Reform the data frame from the matrix
                df <- reshape2::melt(mat, varnames = c("Filename", "GC_Content"), value.name = "Freq")
                df$Filename <- as.character(df$Filename)
                df$Freq <- round(as.numeric(df$Freq), 4)
                df$GC_Content <- as.integer(df$GC_Content)
              }

              key <- unique(df$Filename)
              df$Filename <- labels[df$Filename]
              df$Filename <- factor(df$Filename, levels = unique(df$Filename))
              # Draw the heatmap
              gcPlot <- ggplot(df, aes_string(x = "GC_Content", y = "Filename", fill = "Freq")) +
                geom_tile() +
                scale_x_continuous(expand = c(0, 0)) +
                theme(panel.grid.minor = element_blank(),
                      panel.background = element_blank()) +
                labs(fill = fillLab) +
                scale_fill_gradient2(low = inferno(1, begin = 0.4),
                                     high = inferno(1, begin = 0.9),
                                     midpoint = 0, mid = inferno(1, begin = 0))

              if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

              if(usePlotly){

                gcPlot <- gcPlot +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank())

                t <- getSummary(x)
                t <- t[t$Category == "Per sequence GC content",]
                t$Filename <- factor(labels[t$Filename], levels = levels(df$Filename))
                t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")

                if (missing(pwfCols)) pwfCols <- ngsReports::pwf

                sideBar <- ggplot(t, aes(x = 1, y = Filename, key = key)) +
                  geom_tile(aes_string(fill = "Status")) +
                  scale_fill_manual(values = getColours(pwfCols)) +
                  theme(panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.position="none",
                        axis.title=element_blank(),
                        axis.text=element_blank(),
                        axis.ticks=element_blank())
                sideBar <- suppressMessages(
                  plotly::ggplotly(sideBar, tooltip = c("Filename", "Status"))
                )

                #plot dendrogram
                if(dendrogram && clusterNames){

                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 1)) +
                    scale_x_continuous(expand = c(0,2))

                  gcPlot <- suppressMessages(
                    plotly::subplot(dendro, sideBar, gcPlot, widths = c(0.1,0.1,0.8),
                                    margin = 0, shareY = TRUE) %>%
                      plotly::layout(xaxis3 = list(title = "GC Content (%)"))
                  )
                }
                else{
                  gcPlot <- suppressWarnings(
                    suppressMessages(
                      plotly::subplot(plotly::plotly_empty(), sideBar, gcPlot,
                                      widths = c(0.1,0.1,0.8), margin = 0, shareY = TRUE) %>%
                        plotly::layout(xaxis3 = list(title = "GC Content (%)"),
                                       annotations = list(text = "Filename", showarrow = FALSE,
                                                          textangle = -90))
                    )
                  )
                }

              }
            }
            gcPlot

          }
)
