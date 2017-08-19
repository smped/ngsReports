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
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileName
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
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
plotSummary <- function(x, subset, pwfCols, trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",
                        usePlotly = FALSE, ...){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))
  fillCol <- getColours(pwfCols)

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])

  df <- tryCatch(getSummary(x))

  # Check the pattern contains a capture
  stopifnot(is.logical(trimNames))
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  df$Category <- factor(df$Category, levels = rev(unique(df$Category)))
  df$Status <- factor(df$Status, levels = c("PASS", "WARN", "FAIL"))
  df$StatusNum <- as.integer(df$Status)

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
    sumPlot <- ggplot(df, aes_string(x = "Filename", y = "Category", fill = "StatusNum", key = "Status")) +
      geom_tile(colour = lineCol) +
      geom_vline(xintercept = seq(1.5, nx), colour = lineCol, size = sz) +
      geom_hline(yintercept = seq(1.5, ny), colour = lineCol, size = sz) +
      scale_fill_gradientn(colours = c(fillCol["PASS"], fillCol["WARN"], fillCol["FAIL"]),
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

    plotly::ggplotly(sumPlot, tooltip = c("Filename", "Category", "Status"))

  }
  else{
    sumPlot <- ggplot(df, aes_string(x = "Filename", y = "Category", fill = "Status")) +
      geom_tile(colour = lineCol, size = sz) +
      scale_fill_manual(values = fillCol) +
      labs(x="Filename", y="QC Category") +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    # Add any parameters from dotArgs
    if (!is.null(userTheme)) sumPlot <- sumPlot + userTheme

    sumPlot
  }

}
