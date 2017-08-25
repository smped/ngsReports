#' @title Plot Overrepresented Sequnces As a Heatmap
#'
#' @description Plot Overrepresented Sequences across an entire library
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param pwfCols Object of class \code{PwfCols} to give colours for pass, warning, and fail
#' values in plot
#' @param nSeq \code{numeric}.
#' The number of Sequences to show.
#' Sequences are sorted by \code{max(Percentage)} across all files during selection
#' @param type A regular expression used to subset based on the Possible_Source column
#' @param method Can only take the values \code{"overall"} or \code{"individual"}.
#' Determines whether the top Sequences are selected by the overall ranking (based on Percentage),
#' or whether the top nKmers are selected from each individual file.
#' @param low colour used as the low colour in the heatmap
#' @param high colour used as the high colour in the heatmap
#' @param naCol colour used for missing values
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param clusterNames \code{logical} Cluster heatmap, or not.
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileName
#'
#'
#'
#' @importFrom magrittr %>%
#' @importFrom grDevices rgb
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom ggdendro theme_dendro
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 coord_flip
#'
#' @export
plotOverrepresentedHeatmapPlotly <- function(x,
                                             subset,
                                             pwfCols,
                                             nSeq = 5,
                                             type = ".+",
                                             method = "Individual",
                                             low = rgb(0.2, 0, 0.2), high = rgb(1, 0, 0),
                                             naCol = "grey80",
                                             pattern = "(.+)\\.(fastq|fq).*",
                                             clusterNames = FALSE,
                                             usePlotly = FALSE,
                                             trimNames = FALSE,
                                             dendrogram = FALSE){


   # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.numeric(nSeq))

  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  col <- getColours(pwfCols)

  x <- x[subset]
  df <- tryCatch(Overrepresented_sequences(x))

  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the top Overrepresented Sequences
  nSeq <- as.integer(nSeq)
  if (sum(grepl(type, df$Possible_Source)) == 0){
    warning("No matches to the supplied 'type' could be found. This will be ignored.")
    type <- ".+"
  }

  if (method == "Individual"){
    topSeq <- df %>%
      dplyr::filter(grepl(type, Possible_Source)) %>%
      split(f = .$Filename) %>%
      lapply(dplyr::arrange, desc(Percentage)) %>%
      lapply(dplyr::slice, 1:nSeq) %>%
      dplyr::bind_rows() %>%
      magrittr::extract2("Sequence") %>%
      unique()
  }
  if (method == "Overall"){
    topSeq <- df %>%
      dplyr::filter(grepl(type, Possible_Source)) %>%
      dplyr::group_by(Sequence) %>%
      dplyr::summarise(maxVal = max(Percentage)) %>%
      dplyr::arrange(desc(maxVal)) %>%
      dplyr::slice(1:nSeq) %>%
      magrittr::extract2("Sequence")
  }
  df <- dplyr::filter(df, Sequence %in% topSeq)

  # Create a column which merges the Sequnce & Possible_Source
  df$Seq_Source <- with(df,
                        paste0(Sequence, "\n(", Possible_Source, ")"))

  # Set the Sequence as a factor based on the most prevalent
  seqLevels <- df %>%
    dplyr::arrange(Percentage) %>%
    dplyr::distinct(Seq_Source) %>%
    magrittr::extract2("Seq_Source")
  df$Seq_Source <- factor(df$Seq_Source, levels = seqLevels)

  if (length(seqLevels) < nSeq) {
    message(paste("There is only data in the FASTQC reports for the top",
                  length(seqLevels),"Sequences."))
    nSeq <- length(seqLevels)
  }

  # Quickly add NA values using dcast/melt
  # scale_fill_gradient requires explicit NA values
  df <- df %>%
    reshape2::dcast(Filename~Seq_Source, value.var = "Percentage")

  if(clusterNames){
    xx <- dplyr::select(df, -Filename)
    xx[is.na(xx)] <- 0
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]

  }

  df <- reshape2::melt(df, id.vars = "Filename",
                       variable.name = "Seq_Source", value.name = "Percentage") %>%
    mutate(Filename = factor(Filename, levels = unique(Filename)))

  ORheatmap <- ggplot(data = df, aes(x = Seq_Source, y = Filename, fill = Percentage)) +
    geom_tile() +
    scale_fill_gradient(low = low, high = high, na.value = naCol) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))  +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())

  if(usePlotly){
    t <- dplyr::filter(getSummary(x), Category == "Overrepresented sequences")

    if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
      t <- dplyr::mutate(t, FilenameFull = Filename,
                         Filename = gsub(pattern[1], "\\1", t$Filename),
                         Filename = factor(Filename, levels = unique(Filename)))
    }
    else{
      t <- dplyr::mutate(t,
                         FilenameFull = Filename,
                         Filename = factor(Filename, levels = unique(Filename)))
    }
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
    key <- t$FilenameFull

    sideBar <- ggplot(t, aes(x = 1, y = Filename, key = key)) +
      geom_tile(aes(fill = Status)) +
      scale_fill_manual(values = col) +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position="none",
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    sideBar <- plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))

    ORheatmap <- ORheatmap +
      theme(axis.text = element_blank(), axis.ticks = element_blank())

    #plot dendrogram
    if(dendrogram && clusterNames){
      ggdend <- function(df) {
        ggplot() +
          geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
          theme_dendro()
      }

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 1)) +
        scale_x_continuous(expand = c(0,1))

      dendro <- plotly::ggplotly(dendro) %>%
        plotly::layout(margin = list(b = 0, t = 0))

      ORheatmap <- plotly::subplot(dendro, sideBar, ORheatmap,
                                   widths = c(0.2, 0.1,0.7), margin = 0,
                                   shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "Overrepresented Sequence",
                                     plot_bgcolor = "white"))
    }else{
      ORheatmap <- plotly::subplot(plotly_empty(), sideBar, ORheatmap, widths = c(0.1,0.1,0.8),
                                   margin = 0, shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "Overrepresented Sequence",
                                     plot_bgcolor = "white"),
                       annotations = list(text = "Filename", showarrow = FALSE,
                                          textangle = -90))
    }
  }else{
    ORheatmap <- ORheatmap +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = element_blank()) +
      coord_flip() +
      xlab("Filename") +
      ylab("Sequence")
  }
  ORheatmap
}
