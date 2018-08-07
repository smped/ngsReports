#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' Read totals will be plotted in millions as this is the most common.
#' The raw data from \code{\link{readTotals}} can otherwise be used to manually create a plot.
#'
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
#' @param barCols Colours for duplicated and unique reads.
#' @param expand.x,expand.y Vectors of length 2. Passed to `scale_*_continuous()`
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
#' # Plot the Read Totals without estimated duplicates
#' plotReadTotals(fdl, duplicated = FALSE)
#'
#' @return Returns a ggplot or plotly object
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#'
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
                   labels, millions, barCols, expand.x = c(0, 0), expand.y = c(0, 0), ...){

            df <- readTotals(x)

            stopifnot(is.logical(duplicated))

            # Automatically determine whether to convert to millions
            if (missing(millions)) {
              millions <- ifelse(max(df$Total_Sequences) > 2e06, TRUE, FALSE)
            }
            millions <- millions[1]
            stopifnot(is.logical(millions))
            xlab <- c("Read Totals", "Read Totals (millions)")[millions + 1]

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
            df$Filename <- factor(df$Filename, levels = unique(df$Filename))

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Setup the basic plot in millions or not
            if (millions) df$Total_Sequences <- df$Total_Sequences / 1e06

            # Get the colours for the barplot
            if (missing(barCols)){
              barCols <- c(rgb(0.9, 0.2, 0.2), rgb(0.2, 0.2, 0.8))
            }
            else{
              stopifnot(length(barCols) == (duplicated + 1))
            }

            # Check the axis expansion
            stopifnot(is.numeric(expand.x), is.numeric(expand.y))
            stopifnot(length(expand.x) == 2, length(expand.y) == 2)

            #set left margin
            maxChar <- max(nchar(levels(df$Filename)))
            if(maxChar < 10) l <- 80
            if(maxChar >= 10 & maxChar < 15) l <- 110
            if(maxChar >= 15 & maxChar < 20) l <- 130
            if(maxChar >= 20)  l <- 150

            if (!duplicated){
              df$xmin <- 0
              df$ymax <- as.integer(df[["Filename"]]) + 0.4
              df$ymin <- df[["ymax"]] - 0.9

              rtPlot <- ggplot(df, aes_string(total = "Total_Sequences")) +
                geom_rect(aes_string(xmin = "xmin", xmax = "Total_Sequences",
                                     ymin = "ymin", ymax = "ymax"), fill = barCols[1]) +
                scale_y_continuous(breaks = seq_along(levels(df$Filename)),
                                   labels = levels(df$Filename),
                                   expand = expand.y) +
                scale_x_continuous(expand = expand.x) +
                labs(x = xlab) +
                guides(colour = FALSE) +
                theme_bw()

              if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

              if(usePlotly){

                rtPlot <- suppressMessages(
                  plotly::ggplotly(rtPlot) %>%
                    layout(xaxis = list(title = "Percentage of Total Reads"),
                           margin=list(l=l, b = 50),
                           barmode = "stack")
                )
              }
            }

            if (duplicated){

              stopifnot(bars %in% c("stacked", "adjacent"))

              # Add the information to a joined data.frame
              deDup <- Total_Deduplicated_Percentage(x)
              deDup$Filename <- labels[deDup$Filename]
              deDup$Filename <- factor(deDup$Filename, levels = unique(deDup$Filename))
              colnames(deDup)[which(colnames(deDup) == "Total")] <- "Percentage"


              joinedDf <- dplyr::left_join(deDup, df, by = "Filename")

              if (bars == "adjacent"){


                #Setup the df for plotting
                joinedDf$Deduplicated <- round(joinedDf$Percentage*joinedDf$Total_Sequences/100, 0)
                if (!millions) joinedDf$Deduplicated <- as.integer(round(joinedDf$Deduplicated, 0))

                colnames(joinedDf) <- gsub("Total_Sequences", "Total",colnames(joinedDf))
                joinedDf <- joinedDf[colnames(joinedDf) != "Percentage"]

                joinedDfL <- reshape2::melt(joinedDf, id.vars = "Filename", variable.name = "Type", value.name = c("Total"))

                offsetList <- list(c(0.5, 0.1), c(0.1, -0.4))
                joinedDfL <- split(joinedDfL, joinedDfL$Type)

                joinedDf <- lapply(c(1,2), function(x){

                  df <- joinedDfL[[x]]
                  df$ymin <- as.integer(joinedDf[["Filename"]]) - offsetList[[x]][1]
                  df$ymax <- as.integer(joinedDf[["Filename"]]) - offsetList[[x]][2]
                  df
                })
                joinedDf <- dplyr::bind_rows(joinedDf)
                joinedDf$xmin <- 0
                joinedDf$Type <- factor(joinedDf$Type)

                rtPlot <- ggplot(joinedDf, aes_string(Total = "Total", fill = "Type")) +
                  geom_rect(aes_string(xmin = "xmin", xmax = "Total",
                                       ymin = "ymin", ymax = "ymax")) +
                  scale_fill_manual(values = barCols) +
                  scale_y_continuous(breaks = seq_along(levels(df$Filename)),
                                     labels = levels(df$Filename),
                                     expand = expand.y) +
                  scale_x_continuous(expand = expand.x) +
                  labs(x = xlab) +
                  guides(colour = FALSE) +
                  theme_bw()

                #### need to set colours

              }
              if (bars == "stacked"){

                #Setup the df for plotting
                joinedDf$Unique <- round(joinedDf$Percentage*joinedDf$Total_Sequences/100, 0)


                if (!millions) joinedDf$Unique <- as.integer(round(joinedDf$Unique, 0))
                joinedDf$Duplicated <- joinedDf$Total_Sequences - joinedDf$Unique
                joinedDf <- joinedDf[c("Filename", "Unique", "Duplicated")]
                joinedDf$ymax <- as.integer(joinedDf[["Filename"]]) + 0.4
                joinedDf$ymin <- joinedDf[["ymax"]] - 0.9

                joinedDf <- reshape2::melt(joinedDf, id.vars = c("Filename", "ymin", "ymax"), variable.name = "Type", value.name = "Total")
                joinedDf$Type = factor(joinedDf$Type, levels = c("Duplicated", "Unique"))

                joinedDf <- split(joinedDf, joinedDf$Filename)

                joinedDf <- lapply(joinedDf, function(x){

                  t <- c(1,2)
                  dN <- which(x$Type == "Duplicated")
                  uN <- which(x$Type == "Unique")

                  t[uN] <- 0
                  t[dN] <- x$Total[uN]

                  max <- c(1,2)

                  max[uN] <- x$Total[uN]
                  max[dN] <- x$Total[uN] + x$Total[dN]

                  x$xmin <- t

                  x$xmax <- max
                  x
                })

                joinedDf <- bind_rows(joinedDf)



                rtPlot <- ggplot(joinedDf, aes_string(Total = "Total", fill = "Type", label = "Filename")) +
                  geom_rect(aes_string(xmin = "xmin", xmax = "xmax",
                                       ymin = "ymin", ymax = "ymax")) +
                  scale_fill_manual(values = barCols) +
                  scale_y_continuous(breaks = seq_along(levels(df$Filename)),
                                     labels = levels(df$Filename),
                                     expand = expand.y) +
                  scale_x_continuous(expand = expand.x) +
                  labs(x = xlab) +
                  guides(colour = FALSE) + theme_bw()
              }

              # Add common themes & labels
              if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

              if (usePlotly){

                rtPlot <- suppressMessages(
                  ggplotly(rtPlot, tooltip = c("Total", "fill", "label")) %>%
                    plotly::layout(xaxis = list(title = xlab),
                           margin=list(l=l, b = 50))
                )

              }

            }

            # Draw the plot
            rtPlot

          }
)
