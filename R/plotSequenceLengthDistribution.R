#' @title Plot the Sequence Length Distribution
#'
#' @description Plot the Sequence Length Distribution across one or more FASTQC reports
#'
#' @details
#' This extracts the Sequence Length Distribution from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#' For example, preset axis limits can also be overwritten easily by adding a call to \code{scale_y_continuous}
#' after the call to \code{plotSequenceLengthDistribution}.
#'
#' An alternative interactive plot is available by setting the argument \code{usePlotly = TRUE}.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical}. Output as ggplot2 or plotly object.
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param counts \code{logical} Should distributions be shown as counts or frequencies (percentages)
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} and \code{usePlotly} are \code{FALSE}.
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param ... Used to pass additional attributes to theme()
#' @param expand.x Passed to \code{scale_x_discrete}
#' @param gridlineWidth,gridlineCol Passed to geom_hline and geom_vline to determine
#' width and colour of gridlines
#' @param heatCol The colour scheme for the heatmap
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
#' # Plot as a frequency plot using lines
#' plotSequenceLengthDistribution(fdl)
#'
#' @importFrom dplyr vars
#' @importFrom plotly ggplotly
#' @importFrom plotly layout
#' @importFrom plotly subplot
#' @importFrom viridisLite inferno
#' @import ggplot2
#'
#' @name plotSequenceLengthDistribution
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setGeneric("plotSequenceLengthDistribution",function(x, usePlotly = FALSE, ...){standardGeneric("plotSequenceLengthDistribution")})
#' @aliases plotSequenceLengthDistribution,character
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceLengthDistribution(x, usePlotly,...)
          }
)
#' @aliases plotSequenceLengthDistribution,FastqcFile
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceLengthDistribution(x, usePlotly,...)
          }
)
#' @aliases plotSequenceLengthDistribution,FastqcFileList
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceLengthDistribution(x, usePlotly,...)
          }
)
#' @aliases plotSequenceLengthDistribution,FastqcData
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, ..., expand.x = c(0,0.2)){

            df <- Sequence_Length_Distribution(x)

            if (!length(df)) {
              #stop("No sequence length Module")
              lenPlot <- emptyPlot("No Sequence Length Module Detected")
              
              if(usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
              return(lenPlot)
            }
            
            
            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Add zero counts for lengths either side of the included range
            df <- dplyr::bind_rows(
              lapply(split(df, f = df$Filename),
                     function(x){
                       # Fix global environment error
                       Lower <- NULL
                       dplyr::bind_rows(x,
                                        dplyr::data_frame(Filename = x$Filename[1],
                                                          Lower = c(min(x$Lower) - 1, max(x$Upper) + 1),
                                                          Upper = Lower,
                                                          Length = as.character(Lower),
                                                          Count = 0)
                       )
                     })
            )

            df$Lower <- as.integer(df$Lower)
            df <- reshape2::dcast(df, Filename ~ Lower, value.var = "Count")
            df[is.na(df)] <- 0
            df$Filename <- labels[df$Filename]
            df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Lower", value.name = "Count")
            df$Count <- as.numeric(df$Count)

            # Arrange in position & make the plotting column one called Length
            df <- dplyr::arrange_at(df, vars("Lower"))
            df$Length <- factor(df$Lower, levels = unique(df$Lower))

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            lenPlot <- ggplot(df, aes_string("Length", "Count", colour = "Filename", group = "Filename")) +
              geom_line() +
              facet_wrap(~Filename) +
              labs(y = "Count") +
              scale_x_discrete(expand = expand.x) +
              labs(x = "Sequence Length") +
              theme_bw() +
              theme(legend.position = "none")

            if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

            if(usePlotly){
              lenPlot <- suppressMessages(
                ggplotly(lenPlot, tooltip = c("x", "y"))
              )
              
              lenPlot <- suppressMessages(
                plotly::subplot(plotly::plotly_empty(), lenPlot, widths = c(0.14,0.86)) %>% 
                  layout(xaxis2 = list(title = "Sequence Length (bp)"), yaxis2 = list(title = "Count")))
              
            }
            
            lenPlot

          }
)
#' @aliases plotSequenceLengthDistribution,FastqcDataList
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, counts = FALSE, plotType = "heatmap",
                   cluster = FALSE, dendrogram = FALSE, ...,
                   expand.x = c(0,0.2), gridlineCol = "grey20",
                   gridlineWidth = 0.2, heatCol = inferno(50)){

            df <- Sequence_Length_Distribution(x)

            if (!length(df)) {
              #stop("No sequence length Module")
              lenPlot <- emptyPlot("No Sequence Length Module Detected")
              
              if(usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
              return(lenPlot)
            }
            
            
            # Check for valid plotType
            stopifnot(plotType %in% c("line", "heatmap"))

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                                  names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Add zero counts for lengths either side of the included range
            df <- dplyr::bind_rows(
              lapply(split(df, f = df$Filename),
                     function(x){
                       # Fix global environment error
                       Lower <- NULL
                       dplyr::bind_rows(x,
                                        dplyr::data_frame(Filename = x$Filename[1],
                                                          Lower = c(min(x$Lower) - 1, max(x$Upper) + 1),
                                                          Upper = Lower,
                                                          Length = as.character(Lower),
                                                          Count = 0)
                       )
                     })
            )

            df$Lower <- as.integer(df$Lower)
            df <- reshape2::dcast(df, Filename ~ Lower, value.var = "Count")
            df[is.na(df)] <- 0


            if(cluster){
              xx <- dplyr::select(df, -Filename)
              clus <- as.dendrogram(hclust(dist(xx)))
              row.ord <- order.dendrogram(clus)
              df <- df[row.ord,]
            }
            df$Filename <- labels[df$Filename]
            df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Lower", value.name = "Count")
            df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

            df$Count <- as.numeric(df$Count)
            # Convert the counts to frequencies
            df <- dplyr::bind_rows(
              lapply(split(df, f = df$Filename),
                     function(x){
                       x$Freq <- 100*x$Count / sum(x$Count)
                       x
                     }))

            # Arrange in position & refer to the 'Lower' column as Length
            df <- dplyr::arrange_at(df, vars("Filename", "Lower"))
            df$Length <- factor(df$Lower, levels = unique(df$Lower))

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            if (plotType == "line"){
              if (counts){
                lenPlot <- ggplot(df, aes_string("Length", "Count", colour = "Filename", group = "Filename")) +
                  geom_line() +
                  labs(y = "Count") +
                  theme_bw()
              }
              else{
                lenPlot <- ggplot(df, aes_string("Length", "Freq", colour = "Filename", group = "Filename")) +
                  geom_line() +
                  scale_y_continuous() +
                  labs(y = "Percent (%)") +
                  theme_bw()
              }
            }

            if (plotType == "heatmap"){

              if (counts){
                lenPlot <- ggplot(df, aes_string("Length", "Filename", fill = "Count")) +
                  geom_tile(colour = gridlineCol) +
                  scale_fill_gradientn(colours = heatCol) +
                  scale_y_discrete(labels = labels, expand = c(0, 0))
              }
              else{
                lenPlot <- ggplot(df, aes_string(x = "Length", y ="Filename", fill = "Freq")) +
                  geom_tile(colour = gridlineCol) +
                  labs(fill = "Percent (%)") +
                  scale_fill_gradientn(colours = heatCol) +
                  scale_y_discrete(labels = labels, expand = c(0, 0))
              }
            }

            lenPlot <- lenPlot +
              scale_x_discrete(expand = expand.x) +
              labs(x = "Sequence Length") +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
            if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

            if(usePlotly){
              if(plotType == "line") {
                lenPlot <- lenPlot +
                  theme(legend.position = "none")
                lenPlot <- suppressMessages(
                  ggplotly(lenPlot, tooltip = c("x", "y", "colour"))
                )
              }
              else{

                pwfCols <- ngsReports::pwf
                
                lenPlot <- lenPlot  +
                  theme(axis.ticks.y = element_blank(),
                        axis.text.y = element_blank(),
                        panel.background = element_blank())

                t <- getSummary(x)
                t <- t[t$Category == "Sequence Length Distribution",]
                t$Filename <- labels[t$Filename]
                t$Filename <- factor(t$Filename, levels = levels(df$Filename))
                t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")

                # Set the key for extracting individual files on click
                key <- names(labels[match(levels(df$Filename), labels)])
                sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)

                #plot dendrogram
                if(dendrogram && cluster){

                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0,0.5))
                }
                else{
                  dendro <- plotly::plotly_empty()
                }

                lenPlot <- suppressMessages(
                  plotly::subplot(dendro, sideBar, lenPlot, widths = c(0.1, 0.08, 0.82),
                                  margin = 0.001, shareY = TRUE, titleX = TRUE)
                )
              }

            }

            lenPlot

          }
)
