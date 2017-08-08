#' @title Plot the GC content as a heatmap
#'
#' @description Plot the Per Sequence GC Content for a set of FASTQC files
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
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
#' plotGcHeatmap(fdl)
#'
#' # Using counts
#' plotGcHeatmap(fdl, counts = TRUE)
#'
#' 
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#'
#' @export
plotGcHeatmap <- function(x, subset, counts = FALSE,
                          trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

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
                    Filename = "Observed Mean")

    df <- dplyr::bind_rows(df, mn)
    df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

    gcPlot <- ggplot(df, aes(x = GC_Content, y = Filename, fill = Freq)) +
      labs(fill = "Frequency")

  }
  else{

    # Add the observed mean at each value:
    mn <- dplyr::group_by(df, GC_Content) %>%
      dplyr::summarise(Count = mean(Count)) %>%
      dplyr::mutate(Filename = "Observed Mean")

    df <- dplyr::bind_rows(df, mn)
    df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

    # Initialise the plot using counts
    gcPlot <- ggplot(df, aes(x = GC_Content, y = Filename, fill= Count)) +
      labs(fill = "Count")

  }

  # Add the rest of the plotting detail
  gcPlot <- gcPlot  +
    geom_tile() +
    xlab("GC Content") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0))

  # Draw the plot
  gcPlot

}
