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
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param trimNames \code{logical}. Remove the file suffix from the names displyed in the legend.
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
#' # Draw the default plot
#' plotNContent(fdl)
#'
#'
#'
#' 
#' @importFrom stringr str_detect
#' 
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#'
#' @export
plotNContent <- function(x, subset, pwfCols,
                         trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))
  stopifnot(is.logical(trimNames))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }

  # Get the NContent
  x <- tryCatch(x[subset])
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
  gradCols <- getColours(pwfCols)[1:nCols]
  breaks <- c(0, 5, 20, upr)[1:nCols]

  ggplot(df, aes(x = Base, y = Filename, fill = Percentage)) +
    geom_tile() +
    # scale_fill_gradient2(low = pass, mid = warn, high = fail, midpoint = midpoint) +
    scale_fill_gradientn(colours = gradCols, values = breaks/max(breaks)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}
