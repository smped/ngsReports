#' @title Plot the read totals after de-duplication
#'
#' @description Plot the expected read totals after de-duplication
#'
#' @details
#' This uses the percentages from T\code{Total_Deduplicated_Percentage} and projects
#' them using the existing \code{readTotals}.
#' These numbers are then plotted against the existing readTotals using \code{\link{geom_bar}}
#'
#' Setting \code{bars = "stacked"} gives unique and duplicated reads as a stacked bar graph
#' for each file.
#' Setting \code{bars = "adjacent"} gives the total and deduplicated (i.e. unique) reads as
#' adjacent bars on the graph for each file.
#'
#' Plots can be customised using the standard methods of ggplot2
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical}. If \code{TRUE} the plot will render as an interactive plot in the Viewer pane.
#' Otherwise the plot will be generated using the current plot device.
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param millions \code{logical}. Use Millions of reads as the scale for the y-axis.
#' Unless specified, will be set as TRUE automatically if the highest total is > 2e06.
#' @param bars \code{character} Can take the values \code{bars = "stacked"} or \code{bars = "adjacent"}.
#' Specifying \code{bars = "stacked"} will plot stacked bars with Duplicated and Unique reads.
#' Specifying \code{bars = "adjacent"} will plot Totals before and after de-duplication.
#' @param dupCol The colour to use for the Total values
#' @param uniqCol The colour to use for the Deduplicated values
#' @param ... Used to pass additional attributes to theme()
#'
#' @return Returns a ggplot object, or interactive plotly object
#'
#' @examples
#' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # Draw the plot with stacked bars
#' plotDeduplicatedTotals(fdl)
#'
#' @importFrom grDevices rgb
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom magrittr %>%
#'
#' @export
plotDeduplicatedTotals <- function(x, usePlotly = FALSE, labels,
                                   millions, bars = "stacked",
                                   dupCol = rgb(0.2, 0.2, 0.8), uniqCol = rgb(0.9, 0.2, 0.2),
                                   ...
                                   ){

  rt <- tryCatch(readTotals(x))
  deDup <- tryCatch(Total_Deduplicated_Percentage(x))

  # Automatically determine whether to convert to millions
  if (missing(millions)) {
    millions <- ifelse(max(rt$Total_Sequences) > 2e06, TRUE, FALSE)
  }
  millions <- millions[1]
  stopifnot(is.logical(millions))
  ylab <- c("Read Totals", "Read Totals (millions)")[millions + 1]

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(rt$Filename)),
                        names = unique(rt$Filename))
  }
  else{
    if (!all(unique(rt$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  if (bars == "adjacent"){
    joinedDf <- deDup %>%
      dplyr::rename(Percentage = Total) %>%
      dplyr::left_join(rt, by = "Filename") %>%
      dplyr::mutate(Deduplicated = floor(Percentage*Total_Sequences/100)) %>%
      dplyr::rename(Total = Total_Sequences) %>%
      dplyr::select(Filename, Total, Deduplicated) %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Type", value.name = "Total") %>%
      dplyr::mutate(Total = Total / (10^(millions*6)))

    deDupPlot <- ggplot(joinedDf, aes_string(x = "Filename", y = "Total", fill = "Type")) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = (c(Total = dupCol, Deduplicated = uniqCol))) +
      scale_x_discrete(labels = labels)

  }
  if (bars == "stacked"){

    joinedDf <- deDup %>%
      dplyr::rename(Percentage = Total) %>%
      dplyr::left_join(rt, by = "Filename") %>%
      dplyr::mutate(Unique = floor(Percentage*Total_Sequences/100),
                    Duplicated = Total_Sequences - Unique) %>%
      dplyr::select(Filename, Unique, Duplicated) %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Type", value.name = "Total") %>%
      dplyr::mutate(Total = Total / (10^(millions*6)),
                    Type = factor(Type, levels = c("Duplicated", "Unique")))

    deDupPlot <- ggplot(joinedDf, aes_string(x = "Filename", y = "Total", fill = "Type")) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = (c(Duplicated = dupCol, Unique = uniqCol))) +
      scale_x_discrete(labels = labels)

  }

  # Add the basic layout
  deDupPlot <- deDupPlot +
    labs(y = ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(userTheme)) deDupPlot <- deDupPlot + userTheme

  # Draw the plot
  if (usePlotly) {
    deDupPlot <- plotly::ggplotly(deDupPlot)
  }

  deDupPlot
}
