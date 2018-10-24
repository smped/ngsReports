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
#' \code{plotType = "line"} or \code{plotType = "cumulative"}
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
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Plot as a frequency plot using lines
#' plotSequenceLengthDistribution(fdl)
#' 
#' # Or plot the cumulative value
#' plotSequenceLengthDistribution(fdl, plotType = "cumulative")
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
          function(x, usePlotly = FALSE, labels, plotType = c("line", "cumulative"),
                   ..., expand.x = c(0,0.2)){
              
              df <- Sequence_Length_Distribution(x)
              plotType <- match.arg(plotType)
              
              if (!length(df)) {

                  lenPlot <- emptyPlot("No Sequence Length Module Detected")
                  if(usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
                  return(lenPlot)
              }
              
              labels <- setLabels(df, labels, ...)
              df$Filename <- labels[df$Filename]
              
              # Add zero counts for lengths either side of the included range
              # This is only required if a single value exists
              if (nrow(df) == 1){
                  df <- dplyr::bind_rows(
                      df,
                      dplyr::mutate(df, Lower = Lower - 1, Count = 0),
                      dplyr::mutate(df, Lower = Lower + 1, Count = 0)
                  )
              }
              
              df$Lower <- as.integer(df$Lower)
              df <- dplyr::arrange_at(df, vars("Lower"))
              df <- df[c("Filename", "Length", "Lower", "Count")]
              df$Cumulative <- cumsum(df$Count)
              df$Length <- factor(df$Lower, levels = unique(df$Lower))
              
              # Sort out some plotting parameters
              xLab <- "Sequence Length (bp)"
              yLab <- c(cumulative = "Cumulative Count", line = "Count")[plotType]
              plotY <- c(cumulative = "Cumulative", line = "Count")[plotType]

              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              lenPlot <- ggplot(df, 
                                aes_string("Length", plotY, 
                                           colour = "Filename", 
                                           group = "Filename")) +
                  geom_line() +
                  facet_wrap(~Filename) +
                  labs(x = xLab, y = yLab) +
                  scale_x_discrete(expand = expand.x) +
                  scale_y_continuous(labels = scales::comma) +
                  theme_bw() +
                  theme(legend.position = "none")
              
              if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme
              
              if(usePlotly){
                  lenPlot <- suppressMessages(
                      plotly::ggplotly(lenPlot, tooltip = c("x", "y"))
                  )
                  
                  lenPlot <- suppressMessages(
                      suppressWarnings(
                          plotly::subplot(plotly::plotly_empty(), lenPlot, 
                                          widths = c(0.14,0.86))
                  ))
                  lenPlot <- plotly::layout(lenPlot,
                                            xaxis2 = list(title = xLab), 
                                            yaxis2 = list(title = yLab))
                  
              }
              
              lenPlot
              
          }
)
#' @aliases plotSequenceLengthDistribution,FastqcDataList
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, counts = FALSE, 
                   plotType = c("heatmap", "line", "cumulative"),
                   cluster = FALSE, dendrogram = FALSE, ...,
                   expand.x = c(0,0.2), gridlineCol = "grey20",
                   gridlineWidth = 0.2, heatCol = inferno(50)){
              
              df <- Sequence_Length_Distribution(x)
              
              if (!length(df)) {
                  lenPlot <- emptyPlot("No Sequence Length Module Detected")
                  if(usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
                  return(lenPlot)
              }
              
              # Check for valid plotType
              plotType <- match.arg(plotType)
              
              # Drop the suffix, or check the alternate labels
              labels <- setLabels(df, labels, ...)
              
              # Add zero counts for lengths either side of the included range
              df <- dplyr::bind_rows(
                  lapply(split(df, f = df$Filename),
                         function(x){
                             # Fix global environment error
                             Lower <- NULL
                             dplyr::bind_rows(x,
                                              tibble::tibble(Filename = x$Filename[1],
                                                             Lower = c(min(x$Lower) - 1, max(x$Upper) + 1),
                                                             Upper = Lower,
                                                             Length = as.character(Lower),
                                                             Count = 0)
                             )
                         })
              )
              
              df$Lower <- as.integer(df$Lower)
              df <- dplyr::arrange_at(df, vars("Filename", "Lower"))
              df <- dplyr::group_by(df, Filename)
              df <- dplyr::mutate(df, Cumulative = cumsum(Count))
              if (!counts) {
                  df <- dplyr::mutate(df, 
                                      Cumulative = Cumulative / max(Cumulative),
                                      Freq = Count / sum(Count))
              }
              df <- dplyr::ungroup(df)
              
              # Arrange in position & refer to the 'Lower' column as Length
              df$Length <- factor(df$Lower, levels = unique(df$Lower))
              
              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              if (plotType %in% c("line", "cumulative")){
                  
                  # Decide whether to plot the Count or cumulative sum
                  # and set all labels
                  plotY <- dplyr::case_when(
                      plotType == "cumulative" ~ "Cumulative",
                      plotType == "line" & counts ~ "Count",
                      plotType == "line" & !counts ~ "Freq"
                  )
                  yLab <- dplyr::case_when(
                      plotType == "cumulative" & counts ~ "Cumulative Count",
                      plotType == "cumulative" & !counts ~ "Cumulative (%)",
                      plotType == "line" & counts ~ "Count",
                      plotType == "line" & !counts ~ "Percent (%)"
                  )
                  yLabelFun <- ifelse(counts, scales::comma, scales::percent)
                  
                  df$Filename <- labels[df$Filename]
                  lenPlot <- ggplot(df, aes_string("Length", plotY,
                                                   colour = "Filename",
                                                   group = "Filename")) +
                      geom_line() +
                      labs(y = yLab) +
                      scale_y_continuous(labels = yLabelFun) +
                      theme_bw()
                  
              }
              
              if (plotType == "heatmap"){
                  
                  if (dendrogram && !cluster){
                      message("cluster will be set to TRUE when dendrogram = TRUE")
                      cluster <- TRUE
                  }
                  plotVal <- ifelse(counts, "Count", "Freq")
                  fillLab <- ifelse(counts, "Count", "Percent (%)")
                  fillLabelFun <- ifelse(counts, scales::comma, scales::percent)
                  
                  # Now define the order for a dendrogram if required
                  # This only applies to a heatmap
                  key <- names(labels)
                  if (cluster){
                      clusterDend <- setClusters(df = df[c("Filename", "Lower", plotVal)], 
                                                 rowVal = "Filename", 
                                                 colVal = "Lower", 
                                                 value = plotVal)
                      key <- labels(clusterDend)
                  }
                  # Now set everything as factors
                  df$Filename <- factor(labels[df$Filename], 
                                        levels = labels[key])
                  
                  lenPlot <- ggplot(df, aes_string("Length", "Filename", 
                                                   fill = plotVal)) +
                      geom_tile(colour = gridlineCol) +
                      labs(fill = fillLab) +
                      scale_fill_gradientn(colours = heatCol, labels = fillLabelFun) +
                      scale_y_discrete(labels = labels, expand = c(0, 0))
              }
              
              lenPlot <- lenPlot +
                  scale_x_discrete(expand = expand.x) +
                  labs(x = "Sequence Length") +
                  theme(
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
                      )
              
              if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme
              
              if(usePlotly){
                  if(plotType %in% c("line", "cumulative")) {
                      lenPlot <- lenPlot +
                          theme(legend.position = "none")
                      lenPlot <- suppressMessages(
                          suppressWarnings(
                              plotly::ggplotly(lenPlot, 
                                               tooltip = c("x", "y", "colour"))
                          )
                      )
                  }
                  else{
                      
                      pwfCols <- ngsReports::pwf
                      
                      lenPlot <- lenPlot  +
                          theme(axis.ticks.y = element_blank(),
                                axis.text.y = element_blank(),
                                panel.background = element_blank())
                      
                      status <- getSummary(x)
                      status <- status[status$Category == "Sequence Length Distribution",]
                      status$Filename <- labels[status$Filename]
                      status$Filename <- factor(status$Filename, 
                                                levels = levels(df$Filename))
                      status <- dplyr::right_join(status, unique(df["Filename"]), 
                                                  by = "Filename")
                      # Draw the PWF status as a sideBar
                      sideBar <- makeSidebar(status = status, key = key, 
                                             pwfCols = pwfCols)
                      
                      #plot dendrogram
                      if (dendrogram){
                          dx <- ggdendro::dendro_data(clusterDend)
                          dendro <- ggdend(dx$segments) +
                              coord_flip() +
                              scale_y_reverse(expand = c(0, 0)) +
                              scale_x_continuous(expand = c(0, 0.5))
                          dendro <- plotly::ggplotly(dendro, tooltip = NULL)
                          lenPlot <- lenPlot + ylab("")
                      }
                      else{
                          dendro <- plotly::plotly_empty()
                      }
                      
                      lenPlot <- suppressMessages(
                          suppressWarnings(
                              plotly::subplot(dendro, sideBar, lenPlot, 
                                              widths = c(0.1, 0.08, 0.82),
                                              margin = 0.001, shareY = TRUE, 
                                              titleX = TRUE)
                          )
                      )
                  }
                  
              }
              
              lenPlot
              
          }
)
