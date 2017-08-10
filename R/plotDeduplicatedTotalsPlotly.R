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
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param millions \code{logical}. Use Millions of reads as the scale for the y-axis.
#' Unless specified, will be set as TRUE automatically if the highest total is > 2e06.
#' @param bars \code{character} Can take the values \code{bars = "stacked"} or \code{bars = "adjacent"}.
#' Specifying \code{bars = "stacked"} will plot stacked bars with Duplicated and Unique reads.
#' Specifying \code{bars = "adjacent"} will plot Totals before and after de-duplication.
#' @param col1 The colour to use for the Total values
#' @param col2 The colour to use for the Deduplicated values
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @return Returns a ggplot object.
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
#' # Draw the plot with adjacent bars
#' plotDeduplicatedTotals(fdl)
#'
#' # Draw the plot with stacked bars
#' plotDeduplicatedTotals(fdl, bars = "stacked")
#'
#'
#' @importFrom grDevices rgb
#' @importFrom plotly layout
#' @importFrom plotly subplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme_bw
#'
#' @export
plotDeduplicatedTotalsPlotly <- function(x, subset, millions, bars = "stacked",
                                         col1 = rgb(0.2, 0.2, 0.8), col2 = rgb(0.9, 0.2, 0.2),
                                         trimNames = FALSE, pattern = "(.+)\\.(fastq|fq).*"){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(bars %in% c("adjacent", "stacked"))

  pwfCols <- ngsReports::pwf
  col <- getColours(pwfCols)

  x <- x[subset]
  rt <- tryCatch(readTotals(x, subset = subset, trimNames = trimNames, pattern = pattern))
  deDup <- tryCatch(Total_Deduplicated_Percentage(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    deDup$Filename <- gsub(pattern[1], "\\1", deDup$Filename)
  }

  # Automatically determine whether to convert to millions
  if (missing(millions)) {
    millions <- ifelse(max(rt$Total_Sequences) > 2e06, TRUE, FALSE)
  }
  millions <- millions[1]
  stopifnot(is.logical(millions))
  ylab <- c("Read Totals", "Read Totals (millions)")[millions + 1]

  if (bars == "adjacent"){
    joinedDf <- deDup %>%
      dplyr::rename(Percentage = Total) %>%
      dplyr::left_join(rt, by = "Filename") %>%
      dplyr::mutate(Deduplicated = floor(Percentage*Total_Sequences/100)) %>%
      dplyr::rename(Total = Total_Sequences) %>%
      dplyr::select(Filename, Total, Deduplicated) %>%
      reshape2::melt(id.vars = "Filename", variable.name = "Type", value.name = "Total") %>%
      dplyr::mutate(Total = Total / (10^(millions*6)))

    deDupPlot <- ggplot(joinedDf, aes(x = Filename, y = Total, fill = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = (c(Total = col1, Deduplicated = col2)))

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

    deDupPlot <- ggplot(joinedDf, aes(x = Filename, y = Total, fill = Type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = (c(Duplicated = col1, Unique = col2)))

  }
  key <- unique(deDup["Filename"])
  t <- getSummary(x) %>% dplyr::filter(Category == "Sequence Duplication Levels")
  t <- dplyr::full_join(key["Filename"], t, by = "Filename")
  t$Filename <- with(t, factor(Filename, levels=Filename))

  sideBar <- ggplot(t, aes(x = Filename, y = 1, key = key)) +
    geom_tile(aes(fill = Status)) +
    scale_fill_manual(values = col) +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.title = element_blank())
  sideBar <- plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))

  # Add the basic layout
  deDupPlot <- deDupPlot +
    theme_bw() +
    theme(axis.text.x =element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none")

  # Draw the plot
  plot <- subplot(sideBar, deDupPlot,
                  nrows = 2, shareX = TRUE, heights = c(0.15, 0.85)) %>%
    layout(showlegend = FALSE, yaxis2 = list(title = ylab), margin = list(r = 200))
  plot
}
