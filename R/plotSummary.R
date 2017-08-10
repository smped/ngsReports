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
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param usePlotly \code{logical}. Generate an interactive plot using plotly
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
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#'
#' @export
plotSummary <- function(x, subset, pwfCols, trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",
                        usePlotly = FALSE){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))
  col <- getColours(pwfCols)

  stopifnot(is.logical(trimNames))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])

  df <- tryCatch(getSummary(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  df$Category <- factor(df$Category, levels = rev(unique(df$Category)))
  df$Status <- factor(df$Status, levels = c("PASS", "WARN", "FAIL"))
  df$StatusNum <- as.integer(df$Status)


   if(usePlotly){
     sumPlot <- ggplot(df, aes(x = Filename, y = Category, fill = StatusNum, key = Status)) +
       geom_tile(colour = "black") +
       scale_fill_gradientn(colours = c(col["PASS"], col["WARN"], col["FAIL"]), values = c(0,1)) +
       labs(x="Filename", y="QC Category") +
       scale_x_discrete(expand=c(0,0)) +
       scale_y_discrete(expand=c(0,0)) +
       theme_bw() +
       theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      legend.position = "none")

       plotly::ggplotly(sumPlot, tooltip = c("Filename", "Category", "Status"))

   }else{
     ggplot(df, aes(x = Filename, y = Category, fill = Status)) +
       geom_tile(colour = "black") +
       scale_fill_manual(values = col) +
       labs(x="Filename", y="QC Category") +
       scale_x_discrete(expand=c(0,0)) +
       scale_y_discrete(expand=c(0,0)) +
       theme_bw() +
       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
   }

}
