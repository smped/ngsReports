#' @title Plot the PASS/WARN/FAIL information
#'
#' @description Extract the PASS/WARN/FAIL summaries and plot them
#'
#' @details This uses the standard ggplot2 syntax to create a three colour plot.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly \code{logical}. Generate an interactive plot using plotly
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#'
#' @return A ggplot2 object (\code{usePlotly = FALSE})
#' or an interactive plotly object (\code{usePlotly = TRUE})
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
#' # Check the overall PASS/WARN/FAIL status
#' plotSummary(fdl)
#'
#' # Change theme parameters using the ellipsis
#' library(ggplot2)
#' plotSummary(fdl, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#'
#' # Add a vector of alternative names
#' altNames <- structure(gsub(".fastq", "", fileName(fdl)), names = fileName(fdl))
#' plotSummary(fdl, labels = altNames)
#'
#' # Interactive plot
#' plotSummary(fdl, usePlotly = TRUE)
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom grid unit
#'
#' @export
plotSummary <- function(x, subset, pwfCols, labels, usePlotly = FALSE, ...){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))
  fillCol <- getColours(pwfCols)

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])
  df <- tryCatch(getSummary(x))

  df$Category <- factor(df$Category, levels = rev(unique(df$Category)))
  df$Status <- factor(df$Status, levels = c("PASS", "WARN", "FAIL"))
  df$StatusNum <- as.integer(df$Status)

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  if ("size" %in% names(dotArgs)){
    sz <- dotArgs$size
  }
  else{
    sz <- 0.2
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

  if(usePlotly){
    nx <- length(x)
    ny <- length(unique(df$Category))
    df$Filename <- labels[df$Filename] # Add the new labels
    sumPlot <- ggplot(df, aes_string(x = "Filename", y = "Category", fill = "StatusNum", key = "Status")) +
      geom_tile(colour = lineCol) +
      geom_vline(xintercept = seq(1.5, nx), colour = lineCol, size = sz) +
      geom_hline(yintercept = seq(1.5, ny), colour = lineCol, size = sz) +
      scale_fill_gradientn(colours = c(fillCol["PASS"],
                                       fillCol["WARN"],
                                       fillCol["FAIL"]),
                           values = c(0,1)) +
      labs(x="Filename", y="QC Category") +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(0.01, 0.01, 0.01, 0.04), "npc"),
            legend.position = "none")
    # Add any parameters from dotArgs
    if (!is.null(userTheme)) sumPlot <- sumPlot + userTheme
    suppressMessages(plotly::ggplotly(sumPlot, tooltip = c("Filename", "Category", "Status")))

  }
  else{
    sumPlot <- ggplot(df, aes_string(x = "Filename", y = "Category", fill = "Status")) +
      geom_tile(colour = lineCol, size = sz) +
      scale_fill_manual(values = fillCol) +
      labs(x="Filename", y="QC Category") +
      scale_x_discrete(expand=c(0,0), labels = labels) +
      scale_y_discrete(expand=c(0,0)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    # Add any parameters from dotArgs
    if (!is.null(userTheme)) sumPlot <- sumPlot + userTheme
    sumPlot
  }

}
