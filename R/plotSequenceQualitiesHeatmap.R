#' @title Plot the Per Sequence Quality Scores
#'
#' @description Plot the Per Sequence Quality Scores for a set of FASTQC reports
#'
#' @details Plots the distribution of average sequence quality scores across the set of files.
#' Values can be plotted either as counts (\code{counts = TRUE}) or as frequencies (\code{counts = FALSE}).
#'
#' Any faceting or scale adjustment can be performed after generation of the initial plot,
#' using the standard methods of ggplot2 as desired.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotlyz
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
#' # Draw the defualt plot
#' plotSequenceQualities(fdl)
#'
#' # Get the R1 files
#' r1 <- grepl("R1", fileNames(fdl))
#'
#' # Customise using the R1 subset, plotting counts,
#' # and faceting after the initial function call
#' plotSequenceQualities(fdl, subset = r1, counts = TRUE) +
#'   facet_wrap(~Filename, ncol = 2)
#'
#' 
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' 
#'
#' @export
plotSequenceQualitiesHeatmap <- function(x, subset, counts = FALSE, pwfCols,
                                  trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",
                                  usePlotly = FALSE, clusterNames = FALSE){

  stopifnot(grepl("(Fastqc|character)", class(x)))
  stopifnot(is.logical(trimNames))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  col <- ngsReports::getColours(pwfCols)

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_sequence_quality_scores(x))



  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }
  df <- reshape2::dcast(df, Filename ~ Quality)
  df[is.na(df)] <- 0

  if(clusterNames){
    xx <- dplyr::select(df, -Filename)
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
    df$Filename <- with(df, factor(Filename, levels=Filename))
  }

  df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Quality", value.name = "Count")

  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>%
      dplyr::mutate(Freq = Count / sum(Count)) %>%
      dplyr::ungroup()
    qualPlot <- ggplot(df, aes(x = Quality, y = Filename, fill = Freq))

  }
  else{

    # Initialise the plot using counts
    qualPlot <- ggplot(df, aes(x = Quality, y = Filename, fill = Count))

  }

  qualPlot <- qualPlot + geom_raster() +
    xlab("Mean Sequence Quality Per Read (Phred Score)") +
    scale_fill_gradientn(colours = inferno(150)) +
    ylab("File names") + theme(panel.grid.minor = element_blank(),
                                       panel.background = element_blank())

  if(usePlotly){

    t <- dplyr::filter(getSummary(x), Category == "Per sequence quality scores")
    t <- dplyr::mutate(t, FilenameFull = Filename,
                       Filename = gsub(pattern[1], "\\1", t$Filename),
                       Filename = factor(Filename, levels = unique(df$Filename)))
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
    key <- t$FilenameFull

    d <- ggplot(t, aes(x = 1, y = Filename, key = key, fill = Status)) + geom_tile() +
      scale_fill_manual(values = col) + theme(panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(),
                                                                legend.position="none",
                                                                axis.title=element_blank(),
                                                                axis.text=element_blank(),
                                                                axis.ticks=element_blank())
    d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


    qualPlot <- qualPlot + theme(axis.title.y = element_blank(),
                                          axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank())

    qualPlot <- plotly::subplot(d, qualPlot, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>%
      plotly::layout(xaxis2 = list(title = "Mean Sequence Quality Per Read (Phred Score)"))

    qualPlot
  }else{

    qualPlot

  }

}
