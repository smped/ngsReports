#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' Read totals will be plotted in millions as this is the most common.
#' The raw data from \code{\link{readTotals}} can otherwise be used to manually create a plot.
#'
#' By default, the estimated percentage of dupicated and unique sequences will be drawn.
#' However, this is based on the value shown on FASTQC reports at the top of DeDuplicatedTotals plot,
#' and is known to be inaccurate.
#' As it still gives a good guide as to sequence diversity it is included as the default.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param millions \code{logical}. Use Millions of reads as the scale for the y-axis.
#' Unless specified, will be set as TRUE automatically if the highest total is > 2e06.
#' @param duplicated \code{logical}. Include deduplicated read total estimates to plot charts
#' @param bars If \code{duplicated = TRUE}, show unique and deduplicated reads as "stacked" or "adjacent".
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param dupCol,uniqCol Colours for duplicated and unique reads. Ignored if \code{duplicated = FALSE}
#' @param ... Used to pass additional attributes to theme()
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
#' # Plot the Read Totals showing estimated duplicates
#' plotReadTotals(fdl)
#'
#' # Plot the Read Totals withou estimated duplicates
#' plotReadTotals(fdl, duplicated = FALSE)
#'
#' @return Returns a ggplot or plotly object
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_fill_manual
#'
#' @name plotReadTotals
#' @rdname plotReadTotals-methods
#' @export
setGeneric("plotReadTotals",function(x, usePlotly = FALSE, duplicated = TRUE, bars = "stacked", ...){standardGeneric("plotReadTotals")})
#' @aliases plotReadTotals,character
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "character",
          function(x, usePlotly = FALSE, duplicated = TRUE, bars = "stacked", ...){
            if (length(x) > 1){
              x <- getFastqcData(x)
              plotReadTotals(x, usePlotly, duplicated, bars, ...)
            }
            else{
              stop("plotReadTotals cannot be called on a single FastqcFile")
            }
          }
)
#' @aliases plotReadTotals,FastqcFileList
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, duplicated = TRUE, bars = "stacked", ...){
            x <- getFastqcData(x)
            plotReadTotals(x, usePlotly, duplicated, bars, ...)
          }
)
#' @aliases plotReadTotals,FastqcDataList
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, duplicated = TRUE, bars = "stacked",
                   labels, millions, pwfCols,
                   dupCol = rgb(0.2, 0.2, 0.8), uniqCol = rgb(0.9, 0.2, 0.2),
                   ...){

            df <- tryCatch(readTotals(x))

            stopifnot(is.logical(duplicated))

            # Automatically determine whether to convert to millions
            if (missing(millions)) {
              millions <- ifelse(max(df$Total_Sequences) > 2e06, TRUE, FALSE)
            }
            millions <- millions[1]
            stopifnot(is.logical(millions))
            ylab <- c("Read Totals", "Read Totals (millions)")[millions + 1]

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)),
                                  names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
            df$Filename <- labels[df$Filename]

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Setup the basic plot in millions or not
            if (millions) df$Total_Sequences <- df$Total_Sequences / 1e06


            if (!duplicated){
              # Add the rest of the parameters
              rtPlot <- ggplot(df, aes_string(x = "Filename", y = "Total_Sequences")) +
                geom_bar(stat = "identity") +
                labs(y = ylab) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
              if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

              if(usePlotly){
                rtPlot <- rtPlot +
                  theme(axis.title.x = element_blank())
                rtPlot <- suppressMessages(
                  plotly::ggplotly(rtPlot)
                )
              }
            }

            if (duplicated){

              stopifnot(bars %in% c("stacked", "adjacent"))

              # Add the information to a joined data.frame
              deDup <- tryCatch(Total_Deduplicated_Percentage(x))
              deDup$Filename <- labels[deDup$Filename]
              colnames(deDup)[which(colnames(deDup) == "Total")] <- "Percentage"
              joinedDf <- dplyr::left_join(deDup, df, by = "Filename")

              if (bars == "adjacent"){

                #Setup the df for plotting
                joinedDf$Deduplicated <- joinedDf$Percentage*joinedDf$Total_Sequences/100
                if (!millions) joinedDf$Deduplicated <- as.integer(round(joinedDf$Deduplicated, 0))
                colnames(joinedDf)[which(colnames(joinedDf) == "Total_Sequences")] <- "Total"
                joinedDf <- joinedDf[c("Filename", "Total", "Deduplicated")]
                joinedDf <- reshape2::melt(joinedDf, id.vars = "Filename", variable.name = "Type", value.name = "Total")

                # Draw the plot
                rtPlot <- ggplot(joinedDf, aes_string(x = "Filename", y = "Total", fill = "Type")) +
                  geom_bar(stat = "identity", position = "dodge") +
                  scale_fill_manual(values = (c(Total = dupCol, Deduplicated = uniqCol)))
              }
              if (bars == "stacked"){

                #Setup the df for plotting
                joinedDf$Unique <- joinedDf$Percentage*joinedDf$Total_Sequences/100
                if (!millions) joinedDf$Unique <- as.integer(round(joinedDf$Unique, 0))
                joinedDf$Duplicated <- joinedDf$Total_Sequences - joinedDf$Unique
                joinedDf <- joinedDf[c("Filename", "Unique", "Duplicated")]
                joinedDf <- reshape2::melt(joinedDf, id.vars = "Filename", variable.name = "Type", value.name = "Total")
                joinedDf$Type = factor(joinedDf$Type, levels = c("Duplicated", "Unique"))

                # Draw the plot
                rtPlot <- ggplot(joinedDf, aes_string(x = "Filename", y = "Total", fill = "Type")) +
                  geom_bar(stat = "identity") +
                  scale_fill_manual(values = (c(Duplicated = dupCol, Unique = uniqCol)))
              }

              # Add common themes & labels
              rtPlot <- rtPlot  +
                labs(y = ylab) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
              if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

              if (usePlotly){

                # Sort out the colours
                if (missing(pwfCols)) pwfCols <- ngsReports::pwf
                stopifnot(isValidPwf(pwfCols))

                key <- unique(deDup["Filename"])
                t <- getSummary(x)
                t <- t[t$Category == "Sequence Duplication Levels",]
                t$Filename <- labels[t$Filename]
                t <- dplyr::full_join(key["Filename"], t, by = "Filename")
                t$Filename <- with(t, factor(Filename, levels=Filename))

                sideBar <- ggplot(t, aes(x = Filename, y = 1, key = key)) +
                  geom_tile(aes_string(fill = "Status")) +
                  scale_fill_manual(values = getColours(pwfCols)) +
                  theme(panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.title=element_blank(),
                        axis.text=element_blank(),
                        axis.ticks=element_blank(),
                        legend.title = element_blank())
                sideBar <- suppressMessages(
                  plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))
                )

                # Add the basic layout
                rtPlot <- rtPlot +
                  theme(axis.title.x = element_blank(),
                        axis.text.x =element_blank(),
                        axis.ticks.x=element_blank(),
                        legend.position = "none")
                rtPlot <- suppressMessages(
                  plotly::ggplotly(rtPlot)
                )

                # Make the final layout for the plot
                rtPlot <- plotly::subplot(sideBar, rtPlot,
                                          nrows = 2, shareX = TRUE, heights = c(0.15, 0.85)) %>%
                  plotly::layout(showlegend = FALSE, yaxis2 = list(title = ylab))
              }

            }

            # Draw the plot
            rtPlot

          }
)
