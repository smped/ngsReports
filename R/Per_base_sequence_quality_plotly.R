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
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
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
#' @import tidyr
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#'
#' @export
plotBaseQualitiesPlotly <- function(x, subset, type = "Mean", pwfCols, dendrogram = FALSE, pattern = "(.+)\\.(fastq|fq).*", clusterNames = FALSE, setHeight = "auto"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))
  stopifnot(type %in% c("Mean", "Median"))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])

  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  col <- ngsReports::getColours(pwfCols)


  #load in applicable data

  df <- tryCatch(Per_base_sequence_quality(x))
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start))

  #get longest sequence
  basicStat <- Basic_Statistics(x) %>% dplyr::select(Filename, Longest_sequence)


  #initialize for Mean Base quality
  if(type == "Mean"){
    df <- df %>% dplyr::right_join(basicStat, by = "Filename") %>%
      dplyr::select(Filename, Start, Mean, Longest_sequence)


    #split data into correct lengths and fill NA's
    dfInner <- df %>% split(f = .['Filename']) %>% lapply(function(x){
      dfFill <- data_frame(Start = 1:x$Longest_sequence[1])
      x <- dplyr::right_join(x, dfFill, by = "Start") %>% zoo::na.locf()
    }) %>% dplyr::bind_rows() %>%
      mutate(Start = as.integer(Start)) %>%
      select(-Longest_sequence) %>%
      reshape2::dcast(Filename ~ Start)

    #cluster names true hclust names
    if(clusterNames){
      xx <- dfInner  %>%
        dplyr::select(-Filename)
      xx[is.na(xx)] <- 0
      clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
      row.ord <- order.dendrogram(clus)
      dfInner <- dfInner[row.ord,]
      dfInner$Filename <- with(dfInner, factor(Filename, levels=Filename))
      dfLong <- dfInner %>% tidyr::gather("Start", "Mean", 2:ncol(.))
      dfLong$Mean <- as.integer(dfLong$Mean)
      dfLong$Start <- as.integer(dfLong$Start)

      t <- ngsReports::getSummary(x) %>% dplyr::filter(Category == "Per base sequence quality")
      t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
      t$Filename <- with(t, factor(Filename, levels=Filename))
      key <- t["Filename"]


      p <- ggplot2::ggplot(dfLong, aes(x = Start, y = Filename)) + ggplot2::geom_tile(aes(fill = Mean), color = "white", size = 30) +
        ggplot2::scale_fill_gradientn(colours = c(col["FAIL"], col["FAIL"],col["WARN"], col["WARN"], col["PASS"], col["PASS"]),
                                      values = rescale(c(0,20,20,30,30,40)),
                                      guide = "colorbar", limits=c(0, 40)) +
        ggplot2::theme(panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

      p <- plotly::ggplotly(p) %>% plotly::layout(margin = list(l = 0))

      d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
        ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(),
                                                                  legend.position="none",
                                                                  axis.title=element_blank(),
                                                                  axis.text=element_blank(),
                                                                  axis.ticks=element_blank())
      d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))

      #plot dendrogram
      if(dendrogram){
        ggdend <- function(df) {
          ggplot2::ggplot() +
            ggplot2::geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) + ggdendro::theme_dendro()
        }

        dx <- ggdendro::dendro_data(clus)
        px <- ggdend(dx$segments) + ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand = c(0, 1)) + ggplot2::scale_x_continuous(expand = c(0,1))

        px <- plotly::ggplotly(px) %>% plotly::layout(margin = list(b = 0, t = 0))

        BQheatmap <- plotly::subplot(px, d, p, widths = c(0.2, 0.1,0.7), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis3 = list(title = "Sequencing Cycle"))
      }

      #plot without dendrogram
      if(!dendrogram){
        BQheatmap <- plotly::subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
      }


    }

    #do not cluster Filenames
    if(!clusterNames){
      dfLong <- dfInner %>% tidyr::gather("Start", "Mean", 2:ncol(.))
      dfLong$Mean <- as.integer(dfLong$Mean)
      dfLong$Start <- as.integer(dfLong$Start)

      t <- getSummary(x) %>% dplyr::filter(Category == "Per base sequence quality")
      t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
      t$Filename <- with(t, factor(Filename, levels=Filename))
      key <- t["Filename"]



      p <- ggplot2::ggplot(dfLong, aes(x = Start, y = Filename)) + ggplot2::geom_tile(aes(fill = Mean), color = "white", size = 5) +
        ggplot2::scale_fill_gradientn(colours = c(col["FAIL"], col["FAIL"],col["WARN"], col["WARN"], col["PASS"], col["PASS"]),
                                      values = rescale(c(0,20,20,30,30,40)),
                                      guide = "colorbar", limits=c(0, 40)) +
        ggplot2::theme(panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

      d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
        ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(),
                                                                  legend.position="none",
                                                                  axis.title=element_blank(),
                                                                  axis.text=element_blank(),
                                                                  axis.ticks=element_blank())
      d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


      BQheatmap <- plotly::subplot(d, p, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
    }
  }


  # initialize for Mean Base quality
  if(type == "Median"){
    df <- df %>% dplyr::right_join(basicStat, by = "Filename") %>%
      dplyr::select(Filename, Start, Median, Longest_sequence)


    #split data into correct lengths and fill NA's
    dfInner <- df %>% split(f = .['Filename']) %>% lapply(function(x){
      dfFill <- data_frame(Start = 1:x$Longest_sequence[1])
      x <- dplyr::right_join(x, dfFill, by = "Start") %>% zoo::na.locf()
    }) %>% dplyr::bind_rows() %>%
      mutate(Start = as.integer(Start)) %>%
      select(-Longest_sequence) %>%
      reshape2::dcast(Filename ~ Start)

    #cluster names true hclust names
    if(clusterNames){
      xx <- dfInner %>% dplyr::select(-Filename)
      xx[is.na(xx)] <- 0
      clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
      row.ord <- order.dendrogram(clus)
      dfInner <- dfInner[row.ord,]
      dfInner$Filename <- with(dfInner, factor(Filename, levels=Filename, ordered=TRUE))
      dfLong <- dfInner %>% tidyr::gather("Start", "Median", 2:ncol(.))
      dfLong$Median <- as.integer(dfLong$Median)
      dfLong$Start <- as.integer(dfLong$Start)

      t <- getSummary(x) %>% dplyr::filter(Category == "Per base sequence quality")
      t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
      t$Filename <- with(t, factor(Filename, levels=Filename))
      key <- t["Filename"]



      p <- ggplot2::ggplot(dfLong, aes(x = Start, y = Filename)) + ggplot2::geom_tile(aes(fill = Median), color = "white", size = 30) +
        ggplot2::scale_fill_gradientn(colours = c(col["FAIL"], col["FAIL"],col["WARN"], col["WARN"], col["PASS"], col["PASS"]),
                                      values = rescale(c(0,20,20,30,30,40)),
                                      guide = "colorbar", limits=c(0, 40)) +
        ggplot2::theme(panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

      d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
        ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(),
                                                                  legend.position="none",
                                                                  axis.title=element_blank(),
                                                                  axis.text=element_blank(),
                                                                  axis.ticks=element_blank())
      d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


      #plot dendrogram
      if(dendrogram){
        ggdend <- function(df) {
          ggplot2::ggplot() +
            geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
            ggplot2::labs(x = "", y = "") + ggplot2::theme_minimal() +
            ggplot2::theme(axis.text = element_blank(), axis.ticks = element_blank(),
                           panel.grid = element_blank())
        }

        dx <- ggdendro::dendro_data(clus)
        px <- ggdend(dx$segments) + ggplot2::coord_flip() + ggplot2::scale_y_reverse() + ggplot2::scale_x_continuous(expand = c(0, 1))
        BQheatmap <- plotly::subplot(px, d, p, widths = c(0.2, 0.1,0.7), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis3 = list(title = "Sequencing Cycle"))
      }

      #plot without dendrogram
      if(!dendrogram){
        BQheatmap <- plotly::subplot(d, p, widths = c( 0.1,0.9), margin = 0) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
      }


    }

    #do not cluster Filenames
    if(!clusterNames){
      dfLong <- dfInner %>% tidyr::gather("Start", "Median", 2:ncol(.))
      dfLong$Median <- as.integer(dfLong$Median)
      dfLong$Start <- as.integer(dfLong$Start)



      t <- getSummary(x) %>% dplyr::filter(Category == "Per base sequence quality")
      t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
      t$Filename <- with(t, factor(Filename, levels=Filename))
      key <- t["Filename"]



      p <- ggplot2::ggplot(dfLong, aes(x = Start, y = Filename)) + ggplot2::geom_tile(aes(fill = Median), color = "white", size = 5) +
        ggplot2::scale_fill_gradientn(colours = c(col["FAIL"], col["FAIL"],col["WARN"], col["WARN"], col["PASS"], col["PASS"]),
                                      values = rescale(c(0,20,20,30,30,40)),
                                      guide = "colorbar", limits=c(0, 40)) +
        ggplot2::theme(panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

      d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
        ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(),
                                                                  legend.position="none",
                                                                  axis.title=element_blank(),
                                                                  axis.text=element_blank(),
                                                                  axis.ticks=element_blank())
      d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


      BQheatmap <- plotly::subplot(d, p, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
    }
  }
  #finish plot
  BQheatmap
}

