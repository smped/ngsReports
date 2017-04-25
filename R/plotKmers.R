#' @title Plot Overrepresented Kmers
#'
#' @description Plot Overrepresented Kmers across an entire library
#'
#' @details This plots the top overrepresented Kmers for one or more FASTQC reports.
#' Sequences are coloured in order of first appearance from left to right for easy
#' identification of larger sequences.
#'
#' Returned plots can be easily modified using the standard ggplot2 methods
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param nc \code{numeric}. The number of columns to create in the plot layout
#' @param nKmers \code{numeric}.
#' The number of Kmers to show.
#' Kmers are sorted by \code{max(Obs/Exp_Max)} acros all files to select these
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @return A standard ggplot2 object
#'
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr slice
#' @importFrom dplyr data_frame
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr rename
#' @importFrom magrittr extract2
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @export
plotKmers <- function(x, subset, nc = 2, nKmers = 6,
                      trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(is.numeric(nc))
  stopifnot(is.numeric(nKmers))

  x <- x[subset]
  df <- tryCatch(Kmer_Content(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the top kMers
  nKmers <- as.integer(nKmers)
  topKmers <- dplyr::group_by(df, Sequence) %>%
    dplyr::summarise(maxVal = max(`Obs/Exp_Max`)) %>%
    dplyr::arrange(desc(maxVal)) %>%
    dplyr::slice(1:nKmers) %>%
    magrittr::extract2("Sequence")

  df %<>%
    dplyr::filter(Sequence %in% topKmers) %>%
    dplyr::rename(Base = `Max_Obs/Exp_Position`) %>%
    dplyr::mutate(Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                  Start = as.integer(Start)) %>%
    dplyr::select(Filename, Sequence, `Obs/Exp_Max`, Start)   # The PValue and Count columns are ignored for this plot

  # Set the x-axis for plotting
  # As there will be significant gaps in this data,
  # the bins used need to be obtained from another slot.
  # The most complete will be Per_base_sequence_quality
  # These values can then be incorporated in the final df for accurate plotting & labelling
  refForX <- unique(Per_base_sequence_quality(x)$Base)
  refForX <- dplyr::data_frame(Base = as.character(refForX),
                               Start = gsub("([0-9]*)-[0-9]*", "\\1", Base))
  refForX$Start <- as.integer(refForX$Start)

  # In order to get a line plot, zero values need to be added to the missing positions
  # The above reference scale for X will be used to label the missing values
  zeros <- with(df,
                expand.grid(list(Filename = unique(Filename),
                                 Sequence = unique(Sequence),
                                 `Obs/Exp_Max` = 0,
                                 Start = seq(0, max(Start) + 0.5, by = 0.5)),
                            stringsAsFactors = FALSE))

  # After the bind_rows, duplicate values will exist at some positions
  # Spuriously introduced zeros need to be removed
  df %<>%
    dplyr::bind_rows(zeros) %>%
    dplyr::arrange(Filename, Sequence, Start, desc(`Obs/Exp_Max`)) %>%
    dplyr::distinct(Filename, Sequence, Start, .keep_all = TRUE)

  # Set the Sequence as a factor based on the first position it appears
  # This way the colours will appear in order in the guide as well as the plot
  kMerLevels <- df %>%
    dplyr::filter(`Obs/Exp_Max` != 0) %>%
    dplyr::arrange(Start) %>%
    dplyr::distinct(Sequence) %>%
    magrittr::extract2("Sequence")
  df$Sequence <- factor(df$Sequence, levels = kMerLevels)

  # Now draw the basic plots
  kMerPlot <- ggplot2::ggplot(df,
                              ggplot2::aes(x = Start, y = `Obs/Exp_Max`, colour = Sequence)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~Filename, ncol = nc) +
    ggplot2::scale_x_continuous(breaks = refForX$Start,
                                labels = refForX$Base) +
    ggplot2::theme_bw() +
    ggplot2::ylab(expression(paste(log[2], " Obs/Exp"))) +
    ggplot2::xlab("Position in read (bp)")

  # Check for binned x-axis values to decied whether to rotate x-axis labels
  # This should be clear if there are more than 2 characters in the plotted labels
  binned <- max(nchar(dplyr::filter(refForX, Start %in% df$Start)$Base)) > 2
  if (binned) {
    kMerPlot <- kMerPlot +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  # Draw the plot
  kMerPlot

}
