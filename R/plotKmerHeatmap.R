#' @title Plot Overrepresented Kmers As a Heatmap
#'
#' @description Plot Overrepresented Kmers across an entire library
#'
#' @details Plots the over-represented Kmers as a heatmap.
#' Colours will correspond to the p-values.
#' In the presence of only zero-valued p-values the heatmap will be filled in a single colour,
#' with missing values defined simply as p < 0.
#'
#' If non-zero p-values are present, colours will denote -log10 transformed values.
#' In this case zero values will be transformed using -log10(q) + 1,
#' where q is the minimum of the non-zero p-values.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param nKmers \code{numeric}. The number of Kmers to show.
#' @param method Can only take the values \code{"overall"} or \code{"individual"}.
#' Determines whether the top nKmers are selected by the overall ranking (based on Obs/Exp_Max),
#' or whether the top nKmers are selected from each individual file.
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
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
#' # Not a great example
#' plotKmerHeatmap(fdl)
#'
#' # Find the R1 files with CC barcodes & just plot these
#' ccR1 <- grepl("CC.+R1", fileNames(fdl))
#' plotKmerHeatmap(fdl, subset = ccR1, method = "individual", nKmers = 3)
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_flip
#'
#' @export
plotKmerHeatmap <- function(x, subset, nKmers = 12, method = "overall",
                            pwfCols, naCol = "grey80", flip = TRUE,
                            trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))
  stopifnot(is.logical(trimNames))
  stopifnot(is.numeric(nKmers))
  stopifnot(method %in% c("overall", "individual"))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }

  x <- tryCatch(x[subset])
  df <- tryCatch(Kmer_Content(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the top kMers
  nKmers <- as.integer(nKmers)
  if (method == "individual"){
    topKmers <- df %>%
      split(f = .$Filename) %>%
      lapply(dplyr::arrange, desc(`Obs/Exp_Max`)) %>%
      lapply(dplyr::slice, 1:nKmers) %>%
      dplyr::bind_rows() %>%
      magrittr::extract2("Sequence") %>%
      unique()
  }
  if (method == "overall"){
    topKmers <- df %>%
      dplyr::arrange(desc(`Obs/Exp_Max`)) %>%
      dplyr::distinct(Sequence) %>%
      dplyr::slice(1:nKmers) %>%
      magrittr::extract2("Sequence")
  }
  df <- dplyr::filter(df, Sequence %in% topKmers)

  # Set the Sequence as a factor based on the first position it appears
  # This way the colours will appear in order in the guide as well as the plot
  kMerLevels <- df %>%
    dplyr::arrange(`Max_Obs/Exp_Position`) %>%
    dplyr::distinct(Sequence) %>%
    magrittr::extract2("Sequence")
  df$Sequence <- factor(df$Sequence, levels = rev(kMerLevels))

  if (length(kMerLevels) < nKmers) {
    message(paste("There is only data in the FASTQC reports for the top",
                  length(kMerLevels),"Kmers."))
    nKmers <- length(kMerLevels)
  }

  # Set the p-values to the -log10 scale
  df$PValue <- -log10(df$PValue)
  allInf <- all(is.infinite(df$PValue))
  if (allInf) {
    message(paste("All PValues for the top", nKmers, "Kmers are zero.\n",
                  "This plot will be relatively uninformative."))
    df$PValue <- 0
  }
  else{
    # If there are some finite -log10 transformed p-values,
    # just add 1 to the the infinite ones
    df$PValue[is.infinite(df$PValue)] <- max(df$PValue[is.finite(df$PValue)]) + 1
  }
  # Ensure that all files are included on the plot
  if (trimNames){
    allNames <- gsub(pattern[1], "\\1", fileNames(x))
  }
  else{
    allNames <- fileNames(x)
  }
  df$Filename <- factor(df$Filename, levels = rev(allNames))

  # The inclusion of files without values needs to be worked on
  # Maybe columns should be added during a dcast somewhere
  if (allInf) {
    # First change the missing values to ">0" as this is all the information we have
    # This is done by gaining explicit NA values, then transforming
    heatPlot <- df %>%
      reshape2::dcast(Filename~Sequence, value.var = "PValue") %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Sequence", value.name = "PValue") %>%
      dplyr::mutate(Filename = factor(Filename, levels = rev(allNames))) %>%
      dplyr::mutate(PValue = dplyr::if_else(is.na(PValue), ">0", "=0")) %>%
      ggplot(aes(x =Filename, y = Sequence, fill = PValue)) +
      geom_tile(colour = "grey30", alpha = 0.9) +
      scale_fill_manual(values = c(`=0` = warn, `>0` = naCol)) +
      labs(fill = "PValue")
  }
  else{
    # First get explicit NA values
    heatPlot <- df %>%
      reshape2::dcast(Filename~Sequence, value.var = "PValue") %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Sequence", value.name = "PValue") %>%
      dplyr::mutate(Filename = factor(Filename, levels = rev(allNames))) %>%
      ggplot(aes(x =Filename, y = Sequence, fill = PValue)) +
      geom_tile(colour = "grey30", alpha = 0.9) +
      scale_fill_gradient2(low = getColours(pwfCols)["PASS"],
                                    mid = getColours(pwfCols)["WARN"],
                                    high = getColours(pwfCols)["FAIL"],
                                    na.value = naCol,
                                    midpoint = max(df$PValue)/2) +
      labs(fill = expression(paste(-log[10], "P")))
  }

  heatPlot <- heatPlot +
    theme_bw() +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                     panel.grid = element_blank())

  if (nKmers > length(x) && flip){
    heatPlot <- heatPlot + coord_flip()
  }

  # Draw the plot
  heatPlot

}
