#' @title Plot Overrepresented Kmers As a Heatmap
#'
#' @description Plot Overrepresented Kmers across an entire library
#'
#' @details Plots the over-represented Kmers as a heatmap.
#' Colours will correspond to the p-values.
#' In the presence of only zero-valued p-values the heatmap will be filled in a single colour,
#' with missing values defined simply as p < 0.
#'
#' If non-zero p-values are present, colours will denote -log10 transformed values.
#' In this case zero values will be transformed using -log10(q) + 1,
#' where q is the minimum of the non-zero p-values.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param nKmers \code{numeric}. The number of Kmers to show.
#' @param method Can only take the values \code{"overall"} or \code{"individual"}.
#' Determines whether the top nKmers are selected by the overall ranking (based on Obs/Exp_Max),
#' or whether the top nKmers are selected from each individual file.
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param naCol colour used for missing values
#' #' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
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
#' # Not a great example
#' plotKmerHeatmap(fdl)
#'
#' # Find the R1 files with CC barcodes & just plot these
#' ccR1 <- grepl("CC.+R1", fileName(fdl))
#' plotKmerHeatmap(fdl, subset = ccR1, method = "individual", nKmers = 3)
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_flip
#' @importFrom plotly subplot
#' @importFrom plotly layout
#'
#' @export
plotKmerHeatmap <- function(x, subset, nKmers = 12, method = "overall",
                            pwfCols, labels, naCol = "grey80",
                            usePlotly = FALSE, clusterNames = FALSE,
                            dendrogram = FALSE, ...){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))
  stopifnot(is.numeric(nKmers))
  stopifnot(method %in% c("overall", "individual"))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))
  col <- getColours(pwfCols)


  df <- tryCatch(Kmer_Content(x))

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


  # Get the top kMers
  nKmers <- as.integer(nKmers)
  if (method == "individual"){
    topKmers <- df %>%
      split(f = .$Filename) %>%
      lapply(dplyr::arrange, desc(`Obs/Exp_Max`)) %>%
      lapply(dplyr::slice, 1:nKmers) %>%
      dplyr::bind_rows() %>%
      magrittr::extract2("Sequence") %>%
      unique()
  }
  if (method == "overall"){
    topKmers <- df %>%
      dplyr::arrange(desc(`Obs/Exp_Max`)) %>%
      dplyr::distinct(Sequence) %>%
      dplyr::slice(1:nKmers) %>%
      magrittr::extract2("Sequence")
  }
  df <- dplyr::filter(df, Sequence %in% topKmers)

  # Set the Sequence as a factor based on the first position it appears
  # This way the colours will appear in order in the guide as well as the plot
  kMerLevels <- df %>%
    dplyr::arrange(`Max_Obs/Exp_Position`) %>%
    dplyr::distinct(Sequence) %>%
    magrittr::extract2("Sequence")
  df$Sequence <- factor(df$Sequence, levels = kMerLevels)

  if (length(kMerLevels) < nKmers) {
    message(paste("There is only data in the FASTQC reports for the top",
                  length(kMerLevels),"Kmers."))
    nKmers <- length(kMerLevels)
  }

  # Set the p-values to the -log10 scale
  df$PValue <- -log10(df$PValue)
  allInf <- all(is.infinite(df$PValue))
  if (allInf) {
    message(paste("All PValues for the top", nKmers, "Kmers are zero.\n",
                  "This plot will be relatively uninformative."))
    #if(usePlotly) stop(All kmer P-values equal zero)
    df$PValue <- 0
  }
  else{
    # If there are some finite -log10 transformed p-values,
    # just add 1 to the the infinite ones
    df$PValue[is.infinite(df$PValue)] <- max(df$PValue[is.finite(df$PValue)]) + 1
  }


  df <- reshape2::dcast(df, Filename~Sequence, value.var = "PValue")

  if(clusterNames){
    xx <- dplyr::select(df, -Filename)
    xx[is.na(xx)] <- 0
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
  }
  key <- df$Filename
  df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Sequence", value.name = "PValue")
  df$Filename <- labels[df$Filename]
  df$Filename <- factor(df$Filename, levels = unique(df$Filename))


   # The inclusion of files without values needs to be worked on
  # Maybe columns should be added during a dcast somewhere
  if (allInf) {
    # First change the missing values to ">0" as this is all the information we have
    # This is done by gaining explicit NA values, then transforming

    df <- dplyr::mutate(df, PValue = dplyr::if_else(is.na(PValue), ">0", "=0"))

    heatPlot <-  ggplot(df, aes(x = Sequence, y = Filename, fill = PValue)) +
      geom_tile(colour = "grey30", alpha = 0.9) +
      scale_fill_manual(values = c(`=0` = getColours(pwfCols)["WARN"][[1]], `>0` = naCol)) +
      labs(fill = "PValue")
  }
  else{
    # First get explicit NA values

    heatPlot <- ggplot(df, aes(x = Sequence, y = Filename, fill = PValue)) +
      geom_tile(colour = "grey30", alpha = 0.9) +
      scale_fill_gradient2(low = getColours(pwfCols)["PASS"],
                           mid = getColours(pwfCols)["WARN"],
                           high = getColours(pwfCols)["FAIL"],
                           na.value = naCol,
                           midpoint = max(df$PValue, na.rm = TRUE)/2) +
      labs(fill = expression(paste(-log[10], "P")))
  }

  heatPlot <- heatPlot +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())


  if(usePlotly){
    nx <- length(unique(df$Filename))
    heatPlot <- heatPlot +
      geom_hline(yintercept = seq(1.5, nx), colour = lineCol, size = lineWidth) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank()) +
      labs(fill = "-log(10) P")

    t <- dplyr::filter(getSummary(x), Category == "Kmer Content")
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

    #plot dendrogram
    if(dendrogram && clusterNames){

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 1)) +
        scale_x_continuous(expand = c(0,1))

      heatPlot <- plotly::subplot(dendro, sideBar, heatPlot,
                                  widths = c(0.1,0.1,0.8), margin = 0,
                                  shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "Kmer Sequence",
                                     plot_bgcolor = "white"))
    }else{

      heatPlot <- plotly::subplot(plotly::plotly_empty(),
                                  sideBar,
                                  heatPlot,
                                  widths = c(0.1,0.1,0.8),
                                  margin = 0,
                                  shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "Kmer Sequence"),
                       annotations = list(text = "Filename", showarrow = FALSE,
                                          textangle = -90))
  }

  }else{
    heatPlot <- heatPlot +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  }
  # Draw the plot
  heatPlot
}

