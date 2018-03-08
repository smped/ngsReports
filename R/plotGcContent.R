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
#' @param Fastafile a fasta file contains DNA sequences to generate theoretical GC content
#' @param n number of simulated reads to generate theoretical GC content from \code{Fastafile}
#' @param bp simulated read length to generate theoretical GC content from \code{Fastafile}
#' @param species \code{character} if \code{gcTheory} is \code{TRUE} its must be accompanied by a species
#' Species currently supported can be obtained using \code{mData(gcTheoretical)}
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType Takes values "line" or "heatmap"
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
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
#' # Plot GC content with theoretical GC content generated from a given fasta file
#' plotGcContent(fdl, Fastafile = system.file("extdata","Athaliana.TAIR10.tRNA.fasta",
#' package="ngsReports"))
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
setGeneric("plotGcContent",
           function(x, usePlotly = FALSE, labels,
                    theoreticalGC = TRUE, theoreticalType = "Genome",
                    species = "Hsapiens", GCobject, Fastafile, n = 1e+6, bp = 100, ...){standardGeneric("plotGcContent")})
#' @aliases plotGcContent,character
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "character",
          function(x, usePlotly = FALSE, labels,
                   theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, Fastafile , n = 1e+6, bp = 100, ...){
            x <- getFastqcData(x)
            plotGcContent(x, usePlotly, labels = labels, theoreticalGC = theoreticalGC,
                          theoreticalType = theoreticalType, species = species, GCobject = GCobject, Fastafile = Fastafile, n = n, bp = bp, ...)
          }
)
#' @aliases plotGcContent,FastqcFile
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, labels,
                   theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, Fastafile , n = 1e+6, bp = 100, ...){
            x <- getFastqcData(x)
            plotGcContent(x, usePlotly, labels = labels, theoreticalGC = theoreticalGC,
                          theoreticalType = theoreticalType, species = species, GCobject = GCobject,Fastafile = Fastafile, n = n, bp = bp, ...)
          }
)
#' @aliases plotGcContent,FastqcFileList
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, labels,
                   theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, Fastafile , n = 1e+6, bp = 100, ...){
            x <- getFastqcData(x)
            plotGcContent(x, usePlotly, labels = labels, theoreticalGC = theoreticalGC,
                          theoreticalType = theoreticalType, species = species, GCobject = GCobject, Fastafile = Fastafile, n = n, bp = bp, ...)

          }
)
#' @aliases plotGcContent,FastqcData
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels,
                   theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, Fastafile, n = 1e+6, bp = 100, counts = FALSE, lineCols = c("red", "blue"),
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
              if (!missing(Fastafile)){
                gcTheoryDF <- gcFromFasta(Fastafile,n,bp)
                subTitle <- paste("Theoretical Distribution based on file ",Fastafile)
              } else{
                gcFun <- tryCatch(match.arg(tolower(theoreticalType), c("genomes","transcriptomes")))
                avail <- do.call(gcFun, list(object = GCobject))
                stopifnot(species %in% avail$Name)
                gcTheoryDF <- getGC(GCobject, name = species, type = theoreticalType)  
                names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
                subTitle <- paste("Theoretical Distribution based on the", species, theoreticalType)
              }
              gcTheoryDF$Type <- "Theoretical Distribution"
              gcTheoryDF$Filename <- "Theoretical Distribution"
              gcTheoryDF$Freq <- round(gcTheoryDF$Freq,4)
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
                ggtitle(label = labels, subtitle = c())
              gcPlot <- suppressWarnings(
                suppressMessages(
                  # Try using subplot with plotly_empty to align with the heatmap in the app
                plotly::ggplotly(gcPlot, tooltip = c("GC_Content", value, "Type"))) %>%
                    layout(legend = list(x = 0.622, y = 1))
                )
              gcPlot <- suppressMessages(
                plotly::subplot(plotly::plotly_empty(), gcPlot, widths = c(0.14,0.86)) %>% 
                  layout(xaxis2 = list(title = "GC content (%)"), yaxis2 = list(title = ylab)))
              
            }

            # Draw the plot
            gcPlot
          }
)
#' @aliases plotGcContent,FastqcDataList
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels,
                   theoreticalGC = TRUE, theoreticalType = "Genome",
                   species = "Hsapiens", GCobject, Fastafile, n=1e+6, bp=100, plotType = c("heatmap", "line"), pwfCols,
                   cluster = FALSE, dendrogram = TRUE, ...){

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
              if (!missing(Fastafile)){
                gcTheoryDF <- gcFromFasta(Fastafile,n,bp)
                subTitle <- paste("Theoretical Distribution based on file ",Fastafile)
              } else {
                gcFun <- tryCatch(match.arg(tolower(theoreticalType), c("genomes","transcriptomes","custome")))
                avail <- do.call(gcFun, list(object = GCobject))
                stopifnot(species %in% avail$Name)
                gcTheoryDF <- getGC(GCobject, name = species, type = theoreticalType)  
                names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
                gcTheoryDF$Type <- "Theoretical Distribution"
                gcTheoryDF$Filename <- "Theoretical Distribution"
                subTitle <- paste("Theoretical Distribution based on the", species, theoreticalType)
                gcTheoryDF$Freq <- round(gcTheoryDF$Freq,4)
              }
              gcTheoryDF$Type <- "Theoretical Distribution"
              gcTheoryDF$Filename <- "Theoretical Distribution"
              gcTheoryDF$Freq <- round(gcTheoryDF$Freq,4)
            }
            else{
              gcTheoryDF <- c()
              subTitle <- c()
            }

            # Check for valid plotType arguments
            plotType <- match.arg(plotType)

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
                  x$Freq <- x$Freq - unlist(gcTheoryDF$Freq)
                  x
                }) %>%
                  dplyr::bind_rows()
                fillLab <- "Difference from\nTheoretical GC"
              }
              else{
                fillLab <- "Frequency"
              }

              if(cluster){
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

                sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)


                #plot dendrogram
                if(dendrogram && cluster){

                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0, 0.5))

                }
                else{
                  dendro <- plotly::plotly_empty()
                }

                gcPlot <- suppressMessages(
                  plotly::subplot(dendro, sideBar, gcPlot, widths = c(0.1,0.08,0.82), margin = 0.001, shareY = TRUE) %>%
                    plotly::layout(xaxis3 = list(title = "GC Content (%)"))
                )

              }
            }
            gcPlot

          }
)

#Calculate GC content from a set of DNA sequences provided by a fasta file
gcFromFasta <- function(Fastafile,n=1e+6,bp=100){
  ref <- readDNAStringSet(filepath = Fastafile)
  file.gc = "gc.txt"
  fastqcTheoreticalGC::generateDistn(ref, file = file.gc, n = n, bp = bp)
  gc <- data.frame(read.table(file.gc,header = FALSE,stringsAsFactors = FALSE))
  file.remove(file.gc)
  colnames(gc) <- c("GC_Content","Freq")
  gc$Freq <- gc$Freq/sum(gc$Freq)
  gc
}

