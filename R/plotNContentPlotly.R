#' @title Plot the per base quality as a heatmap
#'
#' @description Plot the Per Base Sequence Quality for a set of FASTQC files
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param type \code{character} Type of quality data to be presented "mean" or "median"
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param pwfcols Object of class \code{\link{Pwfcol}} to give colours for pass, warning, and fail
#' values in plot
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param
#'
#' @return A ggplot2 object
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
#' # The default plot
#' plotGcHeatmap(fdl)
#'
#' # Using counts
#' plotGcHeatmap(fdl, counts = TRUE)
#'
#' @import ggplot2
#' @import scales
#' @import plotly
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#'
#' @export
plotNContentPlotly <- function(x, subset, pwfCols, pattern = "(.+)\\.(fastq|fq).*", clusterNames = TRUE){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- fastqcReports::pwf


  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }

  # Get the NContent
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_base_N_content(x))
  df <- dplyr::rename(df, Percentage = `N-Count`) %>%
    dplyr::mutate(Base = factor(Base, levels = unique(Base)))
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start))

  # Define the colour palette
  col <- getColours(pwfCols)

  basicStat <- Basic_Statistics(x) %>% select(Filename, Longest_sequence)

  df <- df %>% dplyr::right_join(basicStat, by = "Filename") %>%
    dplyr::select(Filename, Start, Percentage, Longest_sequence) %>%
    tidyr::spread(Start, Percentage)

  splitLengths <- df %>% split(f = .['Longest_sequence'])

  dfInner <- splitLengths %>% lapply(function(k){
    k <- k %>% dplyr::select(Filename, dplyr::one_of(as.character(1:k$Longest_sequence[1]))) %>%
      t() %>%
      zoo::na.locf() %>%
      t() %>%
      as.data.frame()
  }) %>% do.call(plyr::rbind.fill, .) %>% magrittr::set_rownames(.$Filename)

  if(clusterNames){
    xx <- dfInner %>% dplyr::select(-Filename)
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    dfInner <- dfInner[row.ord,]
    dfInner$Filename <- with(dfInner, factor(Filename, levels=Filename))
    dfLong <- dfInner %>% gather("Start", "Percentage", 2:ncol(.))
    dfLong$Percentage <- as.numeric(dfLong$Percentage)
    dfLong$Start <- as.numeric(dfLong$Start)

    key <- unique(dfLong["Filename"])
    t <- fastqcReports::getSummary(x) %>% dplyr::filter(Category == "Per base N content")
    t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
    t$Filename <- with(t, factor(Filename, levels=Filename))
    test <- spread(dfLong, Start, Percentage)

    p <- ggplot2::ggplot(dfLong, ggplot2::aes(x = Start, y = Filename, fill = Percentage)) +
      ggplot2::geom_raster() +
      # ggplot2::scale_fill_gradient2(low = pass, mid = warn, high = fail, midpoint = midpoint) +
      ggplot2::scale_fill_gradientn(colours = c(col["PASS"], col["PASS"], col["WARN"], col["WARN"], col["FAIL"], col["FAIL"]),
                                    values = rescale(c(0,5,5,20,20,30)),
                                    guide = "colorbar", limits=c(0, 40), breaks = c(0, 5, 10, 20, 40)) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())

    p <- plotly::ggplotly(p) %>% plotly::layout(margin = list(l = 0))

    d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_raster(aes(fill = Status)) +
      ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(),
                                                                legend.position="none",
                                                                axis.title=element_blank(),
                                                                axis.text=element_blank(),
                                                                axis.ticks=element_blank())
    d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))

  }
  if(!clusterNames){
    dfInner$Filename <- with(dfInner, factor(Filename, levels=Filename))
    dfLong <- dfInner %>% gather("Start", "Percentage", 2:ncol(.))
    dfLong$Percentage <- as.numeric(dfLong$Percentage)
    dfLong$Start <- as.numeric(dfLong$Start)

    key <- unique(dfInner["Filename"])
    t <- fastqcReports::getSummary(x) %>% dplyr::filter(Category == "Per base N content")
    t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
    t$Filename <- with(t, factor(Filename, levels=Filename))

    p <- ggplot2::ggplot(dfLong, ggplot2::aes(x = Start, y = Filename, fill = Percentage)) +
      ggplot2::geom_raster() +
      # ggplot2::scale_fill_gradient2(low = pass, mid = warn, high = fail, midpoint = midpoint) +
      ggplot2::scale_fill_gradientn(colours = c(col["PASS"], col["PASS"], col["WARN"], col["WARN"], col["FAIL"], col["FAIL"]),
                                    values = rescale(c(0,5,5,20,20,30)),
                                    guide = "colorbar", limits=c(0, 40), breaks = c(0, 5, 10, 20, 40)) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())

    p <- plotly::ggplotly(p) %>% plotly::layout(margin = list(l = 0))

    d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_raster(aes(fill = Status)) +
      ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(),
                                                                legend.position="none",
                                                                axis.title=element_blank(),
                                                                axis.text=element_blank(),
                                                                axis.ticks=element_blank())
    d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


  }

  NCheatmap <- plotly::subplot(d, p, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
  NCheatmap
}
