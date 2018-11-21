#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' Read totals will be plotted in millions as this is the most common.
#' The raw data from \code{\link{readTotals}} can otherwise be used to manually
#' create a plot.
#'
#' However, this is based on the value shown on FASTQC reports at the top of 
#' DeDuplicatedTotals plot, and is known to be inaccurate.
#' As it still gives a good guide as to sequence diversity it is included as the
#' default.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param duplicated \code{logical}. Include deduplicated read total estimates 
#' to plot charts
#' @param bars If \code{duplicated = TRUE}, show unique and deduplicated reads 
#' as "stacked" or "adjacent".
#' @param barCols Colours for duplicated and unique reads.
#' @param expand.x Passed to \code{expand_scale(mult = expand.x)} for the x-axis.
#' @param ... Used to pass additional attributes to theme()
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
setGeneric("plotReadTotals",
           function(x, usePlotly = FALSE, labels, duplicated, bars, barCols, 
                    expand.x, ...){
    standardGeneric("plotReadTotals")
})
#' @aliases plotReadTotals,character
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "character",
          function(x, usePlotly = FALSE, labels, duplicated = TRUE, 
                   bars = c("stacked", "adjacent"), barCols = c("red","blue"),
                   expand.x = c(0, 0.02), ...){
              if (length(x) > 1){
                  x <- getFastqcData(x)
                  plotReadTotals(x, usePlotly, labels, duplicated, bars, barCols,
                                 expand.x, ...)
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
          function(x, usePlotly = FALSE, labels, duplicated = TRUE, 
                   bars = c("stacked", "adjacent"), barCols = c("red","blue"), 
                   expand.x = c(0, 0.02), ...){
              x <- getFastqcData(x)
              plotReadTotals(x, usePlotly, labels, duplicated, bars, barCols,
                             expand.x, ...)
          }
)
#' @aliases plotReadTotals,FastqcDataList
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, duplicated = TRUE,
                   bars = c("stacked", "adjacent"), barCols = c("red","blue"),
                   expand.x = c(0, 0.02), ...){
              
              df <- readTotals(x)
              stopifnot(is.logical(duplicated))
              
              # Drop the suffix, or check the alternate labels
              labels <- setLabels(df, labels, ...)
              df$Filename <- labels[df$Filename]
              df$Filename <- factor(df$Filename, levels = unique(df$Filename))
              
              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              # Get the colours for the barplot
              barCols <- tryCatch(barCols[seq_len(duplicated + 1)])
              
              # Check the axis expansion
              stopifnot(is.numeric(expand.x))
              xMax <- max(df$Total_Sequences)
              xLab <- "Read Totals"
              
              if (!duplicated){

                  rtPlot <- ggplot(df, aes_string("Filename", "Total_Sequences")) +
                      geom_bar(stat = "identity", fill = barCols) +
                      scale_y_continuous(labels = scales::comma,
                                         limits = c(0, xMax),
                                         expand = expand_scale(mult = expand.x)) +
                      labs(y = xLab) + 
                      coord_flip() +
                      theme_bw()
                  
                  if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme
                  
                  if(usePlotly) rtPlot <- plotly::ggplotly(rtPlot)

              }
              
              if (duplicated){
                  
                  bars <- match.arg(bars)
                  
                  # Add the information to a joined data.frame
                  deDup <- Total_Deduplicated_Percentage(x)
                  deDup$Filename <- labels[deDup$Filename]
                  deDup$Filename <- factor(deDup$Filename, levels = levels(df$Filename))
                  deDup <- dplyr::rename(deDup, Percentage = Total)
                  df <- dplyr::left_join(deDup, df, by = "Filename")
                  
                  #Setup the df for plotting
                  df$Unique <- df$Percentage*df$Total_Sequences/100
                  df$Unique <- round(df$Unique, 0)
                  df$Duplicated <- df$Total_Sequences - df$Unique
                  df <- df[c("Filename", "Unique", "Duplicated")]
                  df <- tidyr::gather(df, key = "Type", value = "Total",
                                      tidyselect::one_of(c("Unique", "Duplicated")))
                  
                  position <- c(adjacent = "dodge",
                                stacked = "stack")[bars]
                  
                  # The x-axis expansion needs to be reset for this one
                  if (bars == "adjacent") xMax <- max(df$Total)*(1 + expand.x[[1]])
                  
                  # Make the plot
                  rtPlot <- ggplot(df, aes_string("Filename", "Total",
                                                  fill = "Type")) +
                      geom_bar(stat = "identity", position = position) +
                      scale_y_continuous(labels = scales::comma,
                                         limits = c(0, xMax),
                                         expand = expand_scale(mult = expand.x)) +
                      scale_fill_manual(values = barCols) +
                      labs(y = xLab) + 
                      coord_flip() +
                      theme_bw()
                  
                  # Add common themes & labels
                  if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme
                  
                  if (usePlotly){
                      
                      rtPlot <- rtPlot +
                          theme(legend.position = "none")
                      
                      rtPlot <- suppressMessages(
                          suppressWarnings(
                              plotly::ggplotly(rtPlot)
                          )
                      )
                  }
                  
              }
              
              # Draw the plot
              rtPlot
              
          }
)
