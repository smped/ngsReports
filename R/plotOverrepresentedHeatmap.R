#' @title Plot Overrepresented Sequnces As a Heatmap
#'
#' @description Plot Overrepresented Sequences across an entire library
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param nSeq \code{numeric}.
#' The number of Sequences to show.
#' Sequences are sorted by \code{max(Percentage)} across all files during selection
#' @param type A regular expression used to subset based on the Possible_Source column
#' @param method Can only take the values \code{"overall"} or \code{"individual"}.
#' Determines whether the top Sequences are selected by the overall ranking (based on Percentage),
#' or whether the top nKmers are selected from each individual file.
#' @param low colour used as the low colour in the heatmap
#' @param high colour used as the high colour in the heatmap
#' @param naCol colour used for missing values
#' @param flip \code{logical}. Enable a call to \code{coord_flip} to determine the best direction
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @return A standard ggplot2 object
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
#' # Another example which isn't ideal
#' plotOverrepresentedHeatmap(fdl)
#'
#' # Dig a bit deeper
#' plotOverrepresentedHeatmap(fdl, subset = r1, flip = FALSE, nSeq = 10)
#'
#' # Check the top 2 sequences with No Hit from each R1 file
#' plotOverrepresentedHeatmap(fdl, subset = r1, type = "No Hit", nSeq = 2, method = "individual")
#'
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#'
#' @export
plotOverrepresentedHeatmap <- function(x,
                                       subset,
                                       nSeq = 20,
                                       type = ".+",
                                       method = "overall",
                                       low = rgb(0.2, 0, 0.2), high = rgb(1, 0, 0),
                                       naCol = "grey80",
                                       flip = TRUE,
                                       trimNames = TRUE,
                                       pattern = "(.+)\\.(fastq|fq).*"){
  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(is.numeric(nSeq))

  x <- x[subset]
  df <- tryCatch(Overrepresented_sequences(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the top Overrepresented Sequences
  nSeq <- as.integer(nSeq)
  if (sum(grepl(type, df$Possible_Source)) == 0){
    warning("No matches to the supplied 'type' could be found. This will be ignored.")
    type <- ".+"
  }

  if (method == "individual"){
    topSeq <- df %>%
      dplyr::filter(grepl(type, Possible_Source)) %>%
      split(f = .$Filename) %>%
      lapply(dplyr::arrange, desc(Percentage)) %>%
      lapply(dplyr::slice, 1:nSeq) %>%
      dplyr::bind_rows() %>%
      magrittr::extract2("Sequence") %>%
      unique()
  }
  if (method == "overall"){
    topSeq <- df %>%
      dplyr::filter(grepl(type, Possible_Source)) %>%
      dplyr::group_by(Sequence) %>%
      dplyr::summarise(maxVal = max(Percentage)) %>%
      dplyr::arrange(desc(maxVal)) %>%
      dplyr::slice(1:nSeq) %>%
      magrittr::extract2("Sequence")
  }
  df <- dplyr::filter(df, Sequence %in% topSeq)


  # Create a column which merges the Sequnce & Possible_Source
  df$Seq_Source <- with(df,
                        paste0(Sequence, "\n(", Possible_Source, ")"))

  # Set the Sequence as a factor based on the most prevalent
  seqLevels <- df %>%
    dplyr::arrange(Percentage) %>%
    dplyr::distinct(Seq_Source) %>%
    magrittr::extract2("Seq_Source")
  df$Seq_Source <- factor(df$Seq_Source, levels = seqLevels)

  if (length(seqLevels) < nSeq) {
    message(paste("There is only data in the FASTQC reports for the top",
                  length(seqLevels),"Sequences."))
    nSeq <- length(seqLevels)
  }

  # Quickly add NA values using dcast/melt
  # scale_fill_gradient requires explicit NA values
  heatPlot <- df %>%
    reshape2::dcast(Filename~Seq_Source, value.var = "Percentage") %>%
    reshape2::melt(id.vars = "Filename",
                   variable.name = "Seq_Source", value.name = "Percentage") %>%
    ggplot2::ggplot(ggplot2::aes(x = Filename, y = Seq_Source, fill = Percentage)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = low, high = high, na.value = naCol) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::xlab("Filename") +
    ggplot2::ylab("Sequence") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid = ggplot2::element_blank())

  if (nSeq > length(x) && flip){
    heatPlot <- heatPlot +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1))
  }

  # Draw the plot
  heatPlot

}
