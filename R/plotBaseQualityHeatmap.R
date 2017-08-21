#' @title Plot the per base quality as a heatmap
#'
#' @description Plot the Per Base Sequence Quality for a set of FASTQC files
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param type \code{character} Type of quality data to be presented "mean" or "median"
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
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
#' plotBaseQualitiesHeatmap(fdl)
#'
#'
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom dplyr right_join
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#'
#' @export
plotBaseQualitiesHeatmap <- function(x, subset, type = "Mean",
                                    pwfCols, dendrogram = FALSE,
                                    labels, clusterNames = FALSE,
                                    usePlotly = FALSE, ...){
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
    sz <- dotArgs$size
  }
  else{
    sz <- 0.2
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

  #get longest sequence
  basicStat <- Basic_Statistics(x) %>% dplyr::select(Filename, Longest_sequence)

  df <- dplyr::right_join(df, basicStat, by = "Filename")

  #initialize for Mean Base quality
  if(type == "Mean"){
    df <- dplyr::select(df, Filename, Start, Data = Mean, Longest_sequence)
  }else{
    df <- dplyr::select(df, Filename, Start, Data = Median, Longest_sequence)
  }
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
    reshape2::dcast(Filename ~ Start)

  #cluster names true hclust names
  if(clusterNames){
    xx <- dplyr::select(df, -Filename)
    xx[is.na(xx)] <- 0
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
  }


  #melt data back to long format and change filenames to fit labels
  df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Start", value.name = "Data")
  df$Filename <- labels[df$Filename]
  # start heatmap for mean and median
  if(type == "Mean"){
    df <- dplyr::mutate(df, Mean = as.numeric(Data),
                        Start = as.integer(Start),
                        Filename = factor(Filename, levels = unique(Filename)))
    BQheatmap <- ggplot(df, aes(x = Start, y = Filename, fill = Mean)) +
      geom_tile() +
      ngsReports:::scale_fill_pwf(df$Mean, pwfCols, breaks = c(0, 20, 30, 40), passLow = FALSE)
  }else{
    df <- dplyr::mutate(df, Median = as.numeric(Data),
                        Start = as.integer(Start),
                        Filename = factor(Filename, levels = unique(Filename)))
    BQheatmap <- ggplot(df, aes(x = Start, y = Filename, fill = Median)) +
      geom_tile() +
      ngsReports:::scale_fill_pwf(df$Median, pwfCols, breaks = c(0, 20, 30, 40), passLow = FALSE)
  }

  BQheatmap <- BQheatmap +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank())

  if(usePlotly){
    BQheatmap <- BQheatmap + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())

    t <- dplyr::filter(getSummary(x), Category == "Per base sequence quality")
    t$Filename <- labels[t$Filename]
    t <- dplyr::mutate(t, Filename = factor(Filename, levels = unique(df$Filename)))
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
    key <- t$Filename

    sideBar <- ggplot(t, aes(x = 1, y = Filename, key = key)) + geom_tile(aes(fill = Status)) +
      scale_fill_manual(values = col) + theme(panel.grid.minor = element_blank(),
                                              panel.background = element_blank(),
                                              legend.position="none",
                                              axis.title=element_blank(),
                                              axis.text=element_blank(),
                                              axis.ticks=element_blank())
    sideBar <- plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))

    #plot dendrogram
    if(dendrogram){
      ggdend <- function(df) {
        ggplot() +
          geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
          ggdendro::theme_dendro()
      }

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 1)) +
        scale_x_continuous(expand = c(0,1))

      dendro <- plotly::ggplotly(dendro) %>%
        plotly::layout(margin = list(b = 0, t = 0))

      BQheatmap <- plotly::subplot(dendro, sideBar, BQheatmap, widths = c(0.2, 0.1,0.7), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis3 = list(title = "Sequencing Cycle"))
    }
    else{
      BQheatmap <- plotly::subplot(sideBar, BQheatmap, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
    }


  }else{
    BQheatmap
  }
  BQheatmap
}
