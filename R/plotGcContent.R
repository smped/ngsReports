#' @title Plot the Per Sequence GC Content
#'
#' @description Plot the Per Sequence GC Content for a set of FASTQC files
#'
#' @details
#' Draws a plot of all GC Content overlaid across a single plot.
#' The mean from the observed data at each value is overlaid as a dotted black line.
#' Addition of a theoretical distribution has not yet been implemented.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileName
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#'
#' @return A ggplot2 object
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
#' plotGcContent(fdl)
#'
#' # Plot the R1 files using counts
#' r1 <- grepl("R1", fileName(fdl))
#' plotGcContent(fdl, subset = r1 , counts = TRUE)
#'
#'
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom plotly ggplotly
#'
#' @export
plotGcContent <- function(x, subset, counts = FALSE,
                          trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",
                          usePlotly = FALSE){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  x <- x[subset]
  df <- tryCatch(Per_sequence_GC_content(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the correct y-axis label
  ylab <- c("Frequency", "Count")[counts + 1]

  # Remove zero counts
  df <- dplyr::filter(df, Count > 0)

  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>%
      dplyr::mutate(Freq = Count / sum(Count)) %>%
      dplyr::ungroup() %>%
      dplyr::select(Filename, GC_Content, Freq)

    # Add the observed mean at each value:
    mn <- dplyr::group_by(df, GC_Content) %>%
      dplyr::summarise(Freq = mean(Freq)) %>%
      dplyr::mutate(Freq = Freq / sum(Freq), # Normalise to 1 for a distribution
                    Filename = "Mean")

    gcPlot <- ggplot(df, aes(x = GC_Content, y = Freq, colour = Filename)) +
      geom_line() +
      geom_line(data = mn, aes(x = GC_Content, y = Freq),
                         colour = "black", linetype = 2)

  }
  else{

    # Add the observed mean at each value:
    mn <- dplyr::group_by(df, GC_Content) %>%
      dplyr::summarise(Count = mean(Count)) %>%
      dplyr::mutate(Filename = "Mean")

    # Initialise the plot using counts
    gcPlot <- ggplot(df, aes(x = GC_Content, y = Count, colour = Filename)) +
      geom_line() +
      geom_line(data = mn, aes(x = GC_Content, y = Count),
                         colour = "black", linetype = 2)

  }

  # Add the rest of the plotting detail
  gcPlot <- gcPlot +
    ylab(ylab) +
    theme_bw()

  if(usePlotly){
    if(!counts) value <- "Freq"
    else value <- "Count"
    gcPlot <- gcPlot + labs(colour = "")
    gcPlot <- ggplotly(gcPlot, tooltip = c("GC_Content", value))
  }

  # Draw the plot
  gcPlot


}
