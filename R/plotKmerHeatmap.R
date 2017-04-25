#' @title Plot Overrepresented Kmers As a Heatmap
#'
#' @description Plot Overrepresented Kmers across an entire library
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param nKmers \code{numeric}.
#' The number of Kmers to show.
#' Kmers are sorted by \code{max(Obs/Exp_Max)} acros all files to select these
#' @param low colour used as the low colour in the heatmap
#' @param mid colour used as the mid colour in the heatmap
#' @param high colour used as thehigh colour in the heatmap
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
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr if_else
#' @importFrom magrittr extract2
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#'
#' @export
plotKmerHeatmap <- function(x, subset, nKmers = 12,
                            low = rgb(0.2, 0, 0.2), mid = rgb(1, 0, 0), high = rgb(1, 1, 0),
                            trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(is.numeric(nKmers))
  stopifnot(is.numeric(infReplace))

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
  df <- dplyr::filter(df, Sequence %in% topKmers)

  # Set the Sequence as a factor based on the first position it appears
  # This way the colours will appear in order in the guide as well as the plot
  kMerLevels <- df %>%
    dplyr::arrange(`Max_Obs/Exp_Position`) %>%
    dplyr::distinct(Sequence) %>%
    magrittr::extract2("Sequence")
  df$Sequence <- factor(df$Sequence, levels = rev(kMerLevels))

  if (length(kMerLevels) < nKmers) {
    message(paste("There is only data in the FASTQC report for the top",
                  length(kMerLevels),"Kmers."))
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
  df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

  if (allInf) {
    # First change the missing values to ">0" as this is all the information we have
    # This is done by gaining explicit NA values, then transforming
    heatPlot <- df %>%
      reshape2::dcast(Filename~Sequence, value.var = "PValue") %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Sequence", value.name = "PValue") %>%
      dplyr::mutate(PValue = dplyr::if_else(is.na(PValue), ">0", "0")) %>%
      ggplot2::ggplot(ggplot2::aes(x =Filename, y = Sequence, fill = PValue)) +
      ggplot2::geom_tile(colour = "grey30", alpha = 0.9) +
      ggplot2::scale_fill_manual(values = c(`0` = "red", `>0` = "grey50")) +
      ggplot2::labs(fill = "PValue")
  }
  else{
    # First get explicit NA values
    heatPlot <- df %>%
      reshape2::dcast(Filename~Sequence, value.var = "PValue") %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Sequence", value.name = "PValue") %>%
      ggplot2::ggplot(ggplot2::aes(x =Filename, y = Sequence, fill = PValue)) +
      ggplot2::geom_tile(colour = "grey30", alpha = 0.9) +
      ggplot2::scale_fill_gradient2(low = low,
                                    mid = mid,
                                    high = high,
                                    na.value = "grey80",
                                    midpoint = max(df$PValue)/2) +
      ggplot2::labs(fill = expression(paste(-log[10], "PValue")))
  }

  heatPlot <- heatPlot +
    ggplot2::theme_bw() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                     panel.grid = ggplot2::element_blank())

  # Draw the plot
  heatPlot
}
