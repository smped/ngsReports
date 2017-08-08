#' @title Plot Overrepresented Sequnces As a Heatmap
#'
#' @description Plot Overrepresented Sequences across an entire library
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param pwfcols Object of class \code{\link{Pwfcol}} to give colours for pass, warning, and fail
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
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
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
#' # Another example which isn't ideal
#' plotOverrepresentedHeatmap(fdl)
#'
#' # Dig a bit deeper
#' plotOverrepresentedHeatmap(fdl, subset = r1, flip = FALSE, nSeq = 10)
#'
#' # Check the top 2 sequences with No Hit from each R1 file
#' plotOverrepresentedHeatmap(fdl, subset = r1, type = "No Hit", nSeq = 2, method = "individual")
#'
#' 
#' @importFrom stringr str_detect
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#'
#' @export
plotOverrepresentedHeatmapPlotly <- function(x,
                                             subset,
                                             pwfCols,
                                             nSeq = 20,
                                             type = ".+",
                                             method = "Individual",
                                             low = rgb(0.2, 0, 0.2), high = rgb(1, 0, 0),
                                             naCol = "grey80",
                                             pattern = "(.+)\\.(fastq|fq).*",
                                             clusterNames = TRUE){
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
    reshape2::dcast(Filename~Seq_Source, value.var = "Percentage") %>%
    reshape2::melt(id.vars = "Filename",
                   variable.name = "Seq_Source", value.name = "Percentage")
 if(clusterNames){
    dfInner <- df %>% spread(Seq_Source, Percentage)
    xx <- dfInner %>% dplyr::select(-Filename)
    xx[is.na(xx)] <- 0
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    dfInner <- dfInner[row.ord,]
    dfInner$Filename <- with(dfInner, factor(Filename, levels=Filename))
    dfLong <- dfInner %>% tidyr::gather("Seq", "Percentage", 2:ncol(.))
    dfLong$Percentage <- as.numeric(dfLong$Percentage)

    key <- unique(dfLong["Filename"])
    t <- getSummary(x) %>% dplyr::filter(Category == "Overrepresented sequences")
    t <- dplyr::full_join(key, t, by = "Filename")
    t$Filename <- with(t, factor(Filename, levels=Filename))


    p <- dfLong %>%
      ggplot(aes(x = Seq, y = Filename, fill = Percentage)) +
      geom_tile() +
      scale_fill_gradient(low = low, high = high, na.value = naCol) +
      theme_bw() +
      theme(axis.text= element_blank(),
                     axis.ticks= element_blank(),
                     panel.grid = element_blank())

    d <- ggplot(t, aes(x = 1, y = Filename, key = key)) + geom_tile(aes(fill = Status)) +
      scale_fill_manual(values = col) + theme(panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(),
                                                                legend.position="none",
                                                                axis.title=element_blank(),
                                                                axis.text=element_blank(),
                                                                axis.ticks=element_blank())
    d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


    ORSheatmap <- subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Overrepresented sequence"), yaxis = list(title = "Filename"))


 }
  if(!clusterNames){

    key <- unique(df["Filename"])
    t <- getSummary(x) %>% dplyr::filter(Category == "Overrepresented sequences")
    t <- dplyr::full_join(key, t, by = "Filename")
    t$Filename <- with(t, factor(Filename, levels=Filename))


    p <- df %>%
      ggplot(aes(x = Seq_Source, y = Filename, fill = Percentage)) +
      geom_tile() +
      scale_fill_gradient(low = low, high = high, na.value = naCol) +
      theme_bw() +
      theme(axis.text= element_blank(),
                     axis.ticks= element_blank(),
                     panel.grid = element_blank())

    d <- ggplot(t, aes(x = 1, y = Filename, key = key)) + geom_tile(aes(fill = Status)) +
      scale_fill_manual(values = col) + theme(panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(),
                                                                legend.position="none",
                                                                axis.title=element_blank(),
                                                                axis.text=element_blank(),
                                                                axis.ticks=element_blank())
    d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


    ORSheatmap <- subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Overrepresented sequence"), yaxis = list(title = "Filename"))


  }
 # Draw the plot
 ORSheatmap

}

