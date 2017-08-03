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
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#'
#' @export
plotNContentPlotly <- function(x, subset, pwfCols, pattern = "(.+)\\.(fastq|fq).*",
                               clusterNames = FALSE, trimNames = TRUE, dendrogram = FALSE,
                               usePlotly = FALSE){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf


  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }

  # Get the NContent
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_base_N_content(fdl))
  df <- dplyr::rename(df, Percentage = `N-Count`) %>%
    dplyr::mutate(Base = factor(Base, levels = unique(Base)))
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start))

  # Define the colour palette
  col <- getColours(pwfCols)

  basicStat <- Basic_Statistics(x) %>% dplyr::select(Filename, Longest_sequence) %>%
    dplyr::mutate(Filename = gsub(pattern[1], "\\1", .$Filename))
  df <- dplyr::right_join(df, basicStat, by = "Filename")

  #split data into correct lengths and fill NA's
  df <- split(df, f = df['Filename']) %>%
    lapply(function(x){
      dfFill <- data.frame(Start = 1:x$Longest_sequence[1])
      x <- dplyr::right_join(x, dfFill, by = "Start") %>%
        zoo::na.locf()
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Start = as.integer(Start)) %>%
    dplyr::select(-Longest_sequence) %>%
    dplyr::select(Filename, Start, Percentage)


  # Convert from wide to long & set the correct variable types

  #cluster names true hclust names
  if(clusterNames){
    df <- reshape2::dcast(df, Filename ~ Start)
    xx <- dplyr::select(df, -Filename)
    xx[is.na(xx)] <- 0
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
    df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Start", value.name = "Percentage")
  }

    df <- dplyr::mutate(df, Percentage = as.numeric(Percentage),
                        Start = as.integer(Start),
                        Filename = factor(Filename, levels = unique(Filename)))
    Nheatmap <- ggplot2::ggplot(df, ggplot2::aes(x = Start, y = Filename, fill = Percentage))

    Nheatmap <- Nheatmap + ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colours = c(col["PASS"], col["PASS"], col["WARN"], col["WARN"], col["FAIL"], col["FAIL"]),
                                    values = scales::rescale(c(0,5,5,20,20,30)),
                                    guide = "colorbar", limits=c(0, 40), breaks = c(0, 5, 10, 20, 40)) +
      ggplot2::theme(panel.grid.minor = element_blank(),
                     panel.background = element_blank())

  if(usePlotly){
    t <- dplyr::filter(getSummary(x), Category == "Per base N content")
    t <- dplyr::mutate(t, FilenameFull = Filename,
                       Filename = gsub(pattern[1], "\\1", t$Filename),
                       Filename = factor(Filename, levels = unique(df$Filename)))
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
    key <- t$FilenameFull

    sideBar <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
      ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(),
                                                                legend.position="none",
                                                                axis.title=element_blank(),
                                                                axis.text=element_blank(),
                                                                axis.ticks=element_blank())
    sideBar <- plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))

    #plot dendrogram
    if(dendrogram){
      ggdend <- function(df) {
        ggplot2::ggplot() +
          ggplot2::geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) + ggdendro::theme_dendro()
      }

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) + ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand = c(0, 1)) + ggplot2::scale_x_continuous(expand = c(0,1))

      dendro <- plotly::ggplotly(dendro) %>% plotly::layout(margin = list(b = 0, t = 0))

      Nheatmap <- plotly::subplot(dendro, sideBar, Nheatmap, widths = c(0.2, 0.1,0.7), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis3 = list(title = "Sequencing Cycle"))
    }else{
      Nheatmap <- plotly::subplot(sideBar, BQheatmap, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
    }


  }else{
    Nheatmap
  }
  Nheatmap
}
