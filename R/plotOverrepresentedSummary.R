#' @title Plot a summary of OVer-represented Sequences
#'
#' @description Plot a summary of OVer-represented Sequences for a set of FASTQC reports
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param source A regular expression used to group the possible source into those
#' matching the pattern, and "Other".
#' Those matching will be classified using the pattern itself, without any supplied braces.
#' Defaults to \code{source = "(Primer|Adapter)"}
#' @param col1 The colour to use for the Total values
#' @param col2 The colour to use for the Deduplicated values
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileName
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName.
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
#' plotOverrepresentedSummary(fdl)
#'
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom grDevices rgb
#'
#' @export
plotOverrepresentedSummary <- function(x, subset, source = "(Primer|Adapter)",
                                       col1 = rgb(0.2, 0.2, 0.8), col2 = rgb(0.9, 0.2, 0.2),
                                       trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){
  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  x <- x[subset]
  df <- tryCatch(Overrepresented_sequences(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  df <- dplyr::mutate(df,
                      Type = dplyr::if_else(grepl(source, Possible_Source), source, "Other"),
                      Type = gsub("[\\(\\)]", "", Type)) %>%
    dplyr::group_by(Filename, Type) %>%
    dplyr::summarise(Percentage = sum(Percentage))

  ggplot(df, aes(x = Filename, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity") +
    ylab("Overrepresented Sequences (% of Total)") +
    scale_fill_manual(values = c(col1, col2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}
