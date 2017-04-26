#' @title Plot the Per Base N Content
#'
#' @description Plot the Per Base N Content for a set of FASTQC files as a heatmap
#'
#' @details Draws a heatmap based on the N content at every position.
#' Three colours are used (green, yellow, red) set at the values 0, 5 and 20 respectively.
#' According to the FASTQC help page values above 5 will issue a WARN,
#' whilst those above 20 will issue a FAIL
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param pass The colour to use for low percentages
#' @param warn The colour to use at the WARN threshold
#' @param fail The colour to use at the FAIL threshold
#' @param trimNames \code{logical}. Remove the file suffix from the names displyed in the legend.
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
#' fileList <- system.file("extdata", fileList, package = "fastqcReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # Draw the default plot
#' plotNContent(fdl)
#'
#'
#'
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#'
#' @export
plotNContent <- function(x, subset,
                         pass = rgb(0, 0.8, 0), warn = rgb(0.9, 0.9, 0.2),
                         fail = rgb(0.8, 0.2, 0.2),
                         trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  # Get the NContent
  x <- x[subset]
  df <- tryCatch(Per_base_N_content(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  df <- dplyr::rename(df, Percentage = `N-Count`) %>%
    dplyr::mutate(Base = factor(Base, levels = unique(Base)))

  # Define the colour palette
  upr <- max(df$Percentage)
  nCols <- findInterval(upr, c(0, 5, 20, Inf)) + 1
  gradCols <- c(pass, warn, fail, rgb(1,1,1))[1:nCols]
  breaks <- c(0, 5, 20, upr)[1:nCols]

  ggplot2::ggplot(df, ggplot2::aes(x = Base, y = Filename, fill = Percentage)) +
    ggplot2::geom_tile() +
    # ggplot2::scale_fill_gradient2(low = pass, mid = warn, high = fail, midpoint = midpoint) +
    ggplot2::scale_fill_gradientn(colours = gradCols, values = breaks/max(breaks)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

}
