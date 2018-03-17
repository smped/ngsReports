#' @title Plot the Per Sequence Quality Scores
#'
#' @description Plot the Per Sequence Quality Scores for a set of FASTQC reports
#'
#' @details Plots the distribution of average sequence quality scores across the set of files.
#' Values can be plotted either as counts (\code{counts = TRUE}) or as frequencies (\code{counts = FALSE}).
#'
#' Any faceting or scale adjustment can be performed after generation of the initial plot,
#' using the standard methods of ggplot2 as desired.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param warn,fail The default values for warn and fail are 5 and 10 respectively (i.e. precentages)
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param alpha set alpha for line graph bounds
#' @param lineWidth,lineCol Passed to geom_hline and geom_vline to determine
#' width and colour of gridlines
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#'
#' @return A standard ggplot2 object, or an interactive plotly object
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
#' # The default plot
#' plotSequenceQualities(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fileName(fdl))
#' plotSequenceQualities(fdl[r1])
#'
#'
#' @importFrom dplyr vars funs
#' @importFrom stats hclust dist
#' @import ggplot2
#'
#' @name plotSequenceQualities
#' @rdname plotSequenceQualities-methods
#' @export
setGeneric("plotSequenceQualities",function(x, usePlotly = FALSE,
                                         ...){standardGeneric("plotSequenceQualities")})
#' @aliases plotSequenceQualities,character
#' @rdname plotSequenceQualities-methods
#' @export
setMethod("plotSequenceQualities", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceQualities(x, usePlotly,...)
          }
)
#' @aliases plotSequenceQualities,FastqcFile
#' @rdname plotSequenceQualities-methods
#' @export
setMethod("plotSequenceQualities", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceQualities(x, usePlotly,...)
          }
)
#' @aliases plotSequenceQualities,FastqcFileList
#' @rdname plotSequenceQualities-methods
#' @export
setMethod("plotSequenceQualities", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceQualities(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceQualities,FastqcData
#' @rdname plotSequenceQualities-methods
#' @export
setMethod("plotSequenceQualities", signature = "FastqcData",
          function(x, usePlotly = FALSE, counts = FALSE, warn = 30, fail = 20, pwfCols,
                   alpha = 0.2, ...){

            df <- Per_sequence_quality_scores(x)

            if(length(df)){
          
            # Sort out the colours
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))

            minQ <- min(df$Quality)

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # make Ranges for rectangles and set alpha
            pwfCols <- setAlpha(pwfCols, alpha)
            rects <- data_frame(ymin = 0,
                                ymax = max(df$Count),
                                xmin = c(0, fail, warn),
                                xmax = c(fail, warn, 41),
                                Status = c("FAIL", "WARN", "PASS"))


            if (!counts){
              Count <- NULL # To avoid NOTE messages in R CMD check
              # Summarise to frequencies & initialise the plot
              df <- dplyr::group_by(df, Filename)
              df <- dplyr::mutate(df, Frequency = Count / sum(Count))
              df <- dplyr::ungroup(df)
              df$Frequency <- round(df$Frequency, 3)
              rects$ymax <- max(df$Frequency)

              qualPlot <- ggplot(df) +
                geom_rect(data = rects,
                          aes_string(xmin = "xmin", xmax = "xmax",
                                     ymin = "ymin", ymax = "ymax", fill = "Status")) +
                geom_line(aes_string(x = "Quality", y = "Frequency", colour = "Filename"))
              
              yax <- "Frequency"


            }
            else{
              # Initialise the plot using counts
              qualPlot <- ggplot(df) +
                geom_rect(data = rects,
                          aes_string(xmin = "xmin", xmax = "xmax",
                                     ymin = "ymin", ymax = "ymax", fill = "Status")) +
                geom_line(aes_string(x = "Quality", y = "Count", colour = "Filename"))
              
              yax <- "Count"

            }

            qualPlot <- qualPlot +
              scale_fill_manual(values = getColours(pwfCols))  +
              scale_y_continuous(limits = c(0, rects$ymax[1]), expand =c(0, 0)) +
              scale_x_continuous(expand = c(0, 0)) +
              scale_colour_discrete()  +
              facet_wrap(~Filename) +
              labs(x = "Mean Sequence Quality Per Read (Phred Score)") +
              guides(fill = FALSE) +
              theme_bw()

            if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

            if(usePlotly){
              if (counts) ymax <- max(df$Count)

              qualPlot <- qualPlot + theme(legend.position = "none")
              qualPlot <- suppressMessages(
                plotly::ggplotly(qualPlot,
                                 hoverinfo = c("x", "y", "colour"))
              )
              
              qualPlot <- suppressMessages(
                plotly::subplot(plotly::plotly_empty(), qualPlot, widths = c(0.14,0.86)) %>% 
                  layout(xaxis2 = list(title = "Mean Sequence Quality Per Read (Phred Score)"), yaxis2 = list(title = yax)))
              
              
              # Set the hoverinfo for bg rectangles to the vertices only,
              # This will effectively hide them
              qualPlot$x$data[[1]]$hoveron <- "points"
              qualPlot$x$data[[2]]$hoveron <- "points"
              qualPlot$x$data[[3]]$hoveron <- "points"
            }
            
            
            }
            
            else{
              qualPlot <- emptyPlot("Per Base N Content Module is missing from the input")
              if(usePlotly) qualPlot <- ggplotly(nPlot, tooltip = "")
            }
            

            # Draw the plot
            qualPlot

          }

)
#' @aliases plotSequenceQualities,FastqcDataList
#' @rdname plotSequenceQualities-methods
#' @export
setMethod("plotSequenceQualities", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, counts = FALSE, pwfCols,
                   labels, plotType = "heatmap", dendrogram = FALSE,
                   cluster = FALSE, lineCol = "grey20",
                   lineWidth = 0.2, alpha = 0.1, warn = 30, fail = 20, ...){

            # Read in data
            df <- Per_sequence_quality_scores(x)
            
            if(length(df)){

            # Sort out the colours
            if(base::missing(pwfCols)) pwfCols <- ngsReports::pwf

            # Drop the suffix, or check the alternate labels
            if(base::missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)),
                                  names = unique(df$Filename))
            }else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            if(plotType == "heatmap"){
              df <- reshape2::dcast(df, Filename ~ Quality, value.var = "Count")
              df[is.na(df)] <- 0

              if(cluster){
                xx <- dplyr::select(df, -Filename)
                clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
                row.ord <- order.dendrogram(clus)
                df <- df[row.ord,]
              }

              key <- df$Filename
              df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Quality", value.name = "Count")
              df$Filename <- labels[df$Filename]
              df$Filename <- factor(df$Filename, levels = unique(df$Filename))

              if (!counts){
                Count <- NULL # To avoid NOTE messages in R CMD check
                # Summarise to frequencies & initialise the plot
                df <- dplyr::group_by(df, Filename)
                df <- dplyr::mutate(df, Frequency = Count / sum(Count))
                df <- dplyr::ungroup(df)
                df$Frequency <- round(df$Frequency, 3)
                qualPlot <- ggplot(df, aes_string(x = "Quality", y = "Filename", fill = "Frequency"))
              }
              else{
                qualPlot <- ggplot(df, aes_string(x = "Quality", y = "Filename", fill = "Count"))
              }

              qualPlot <- qualPlot +
                geom_tile(colour = lineCol, size = lineWidth) +
                xlab("Mean Sequence Quality Per Read (Phred Score)") +
                scale_fill_gradientn(colours = inferno(150)) +
                ylab("File names") +
                theme(panel.grid.minor = element_blank(),
                      panel.background = element_blank())

              if(usePlotly){
                # Add lines and remove axis data
                qualPlot <- qualPlot +
                  theme(axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank())

                t <- getSummary(x)
                t <- t[t$Category == "Per sequence quality scores",]
                t$Filename <- labels[t$Filename]
                t <- dplyr::mutate(t, Filename = factor(Filename, levels = levels(df$Filename)))
                t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")

                # Make sidebar
                sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)

                #plot dendrogram
                if(dendrogram && cluster){

                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0,0.5))

                  qualPlot <- suppressMessages(
                    plotly::subplot(dendro, sideBar, qualPlot,
                                    widths = c(0.1,0.08,0.82), margin = 0.001,
                                    shareY = TRUE)
                  )
                  qualPlot <- suppressMessages(
                    plotly::layout(qualPlot,
                                   xaxis3 = list(title = "Mean Sequence Quality Per Read (Phred Score)",
                                                 plot_bgcolor = "white"))
                  )
                }
                else{

                  qualPlot <- suppressMessages(
                    plotly::subplot(plotly::plotly_empty(),
                                    sideBar,
                                    qualPlot,
                                    widths = c(0.1,0.08,0.82),
                                    margin = 0.01,
                                    shareY = TRUE)
                  )
                  qualPlot <- suppressMessages(
                    plotly::layout(qualPlot,
                                   xaxis3 = list(title = "Mean Sequence Quality Per Read (Phred Score)"),
                                   annotations = list(text = "Filename", showarrow = FALSE,
                                         textangle = -90))
                  )
                }
              }
              qualPlot
            }
            if(plotType == "line"){
              # make Ranges for rectangles and set alpha
              pwfCols <- setAlpha(pwfCols, alpha)
              rects <- data_frame(ymin = 0,
                                  ymax = 0,
                                  xmin = c(0, fail, warn),
                                  xmax = c(fail, warn, 41),
                                  Status = c("FAIL", "WARN", "PASS"))


              if (!counts){
                Count <- NULL # To avoid NOTE messages in R CMD check
                # Summarise to frequencies & initialise the plot
                df <- dplyr::group_by(df, Filename)
                df <- dplyr::mutate(df, Frequency = Count / sum(Count))
                df <- dplyr::ungroup(df)
                rects$ymax <- max(df$Frequency)
                qualPlot <- ggplot(df) +
                  geom_rect(data = rects,
                            aes_string(xmin = "xmin", xmax = "xmax",
                                       ymin = "ymin", ymax = "ymax", fill = "Status")) +
                  geom_line(aes_string(x = "Quality", y = "Frequency", colour = "Filename"))

              }
              else{
                # Initialise the plot using counts
                qualPlot <- ggplot(df) +
                  geom_rect(data = rects,
                            aes_string(xmin = "xmin", xmax = "xmax",
                                       ymin = "ymin", ymax = "ymax", fill = "Status")) +
                  geom_line(aes_string(x = "Quality", y = "Count", colour = "Filename"))

              }

              qualPlot <- qualPlot +
                scale_fill_manual(values = getColours(pwfCols))  +
                scale_y_continuous(limits = c(0, rects$ymax[1]), expand =c(0, 0)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_colour_discrete() +
                labs(x = "Mean Sequence Quality Per Read (Phred Score)") +
                guides(fill = FALSE) +
                theme_bw()

              if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

              if(usePlotly){
                if (counts) ymax <- max(df$Count)

                qualPlot <- qualPlot + theme(legend.position = "none")
                qualPlot <- suppressMessages(
                  plotly::ggplotly(qualPlot,
                                   hoverinfo = c("x", "y", "colour"))
                )
                # Set the hoverinfo for bg rectangles to the vertices only,
                # This will effectively hide them
                qualPlot$x$data[[1]]$hoveron <- "points"
                qualPlot$x$data[[2]]$hoveron <- "points"
                qualPlot$x$data[[3]]$hoveron <- "points"
              }}
            
            }
            else{
              qualPlot <- emptyPlot("Per Base N Content Module is missing from the input")
              if(usePlotly) qualPlot <- ggplotly(qualPlot, tooltip = "")
            }
            
            qualPlot
            }

)
