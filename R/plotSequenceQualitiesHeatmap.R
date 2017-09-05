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
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#'
#' @return A ggplot2 object
#'
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom viridisLite inferno
#'
#' @export
plotSequenceQualitiesHeatmap <- function(x, subset, labels, counts = FALSE, pwfCols,
                                         usePlotly = FALSE, dendrogram = FALSE,
                                         clusterNames = FALSE, ...){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  col <- ngsReports::getColours(pwfCols)

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_sequence_quality_scores(x))

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  if ("size" %in% names(dotArgs)){
    lineWidth <- dotArgs$size
  }
  else{
    lineWidth <- 0.2
  }
  if ("colour" %in% names(dotArgs) || "color" %in% names(dotArgs)){
    i <- which(names(dotArgs) %in% c("colour", "color"))
    lineCol <- dotArgs[[i]]
  }
  else{
    lineCol <- "grey20"
  }
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])


  df <- reshape2::dcast(df, Filename ~ Quality)
  df[is.na(df)] <- 0

  if(clusterNames){
    xx <- dplyr::select(df, -Filename)
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
  }

  key <- df$Filename
  df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Quality", value.name = "Count")
  df$Filename <- labels[df$Filename]

  if (!counts){

    # Summarise to frequencies
    df <- dplyr::mutate(df, Value = Count / sum(Count),
                    Filename = factor(Filename, levels = unique(Filename)))


  }else{
    df <- mutate(df, Value = Count,
                 Filename = factor(Filename, levels = unique(Filename)))
  }

  qualPlot <- ggplot(df, aes(x = Quality, y = Filename, fill = Value)) +
    geom_tile(colour = lineCol, size = lineWidth) +
    xlab("Mean Sequence Quality Per Read (Phred Score)") +
    scale_fill_gradientn(colours = inferno(150)) +
    ylab("File names") +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank())

  if(usePlotly){
    nx <- length(x)
    qualPlot <- qualPlot +
      geom_hline(yintercept = seq(1.5, nx), colour = lineCol, size = lineWidth)

    t <- dplyr::filter(getSummary(x), Category == "Per sequence quality scores")
    t$Filename <- labels[t$Filename]
    t <- dplyr::mutate(t, Filename = factor(Filename, levels = unique(df$Filename)))
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")


    sideBar <- ggplot(t, aes(x = 1, y = Filename, key = key, fill = Status)) +
      geom_tile() +
      geom_hline(yintercept = seq(1.5, nx), colour = lineCol, size = lineWidth) +
      scale_fill_manual(values = col) +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position="none",
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())

    sideBar <- plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))

    qualPlot <- qualPlot +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

    #plot dendrogram
    if(dendrogram && clusterNames){

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 1)) +
        scale_x_continuous(expand = c(0,1))

      qualPlot <- plotly::subplot(dendro, sideBar, qualPlot,
                                  widths = c(0.1,0.1,0.8), margin = 0,
                                  shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "Mean Sequence Quality Per Read (Phred Score)",
                                     plot_bgcolor = "white"))
    }else{

      qualPlot <- plotly::subplot(plotly::plotly_empty(),
                                  sideBar,
                                  qualPlot,
                                  widths = c(0.1,0.1,0.8),
                                  margin = 0,
                                  shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "Mean Sequence Quality Per Read (Phred Score)"),
                       annotations = list(text = "Filename", showarrow = FALSE,
                                          textangle = -90))
    }
  }
  qualPlot

}
