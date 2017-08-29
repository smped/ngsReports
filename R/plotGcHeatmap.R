#' @title Plot GC Content for each file
#'
#' @description Plot the GC content for each file
#'
#' @details makes interactive heatmap for \code{GC content}
#'
#' For large datasets, subsetting by R1 or R2 reads may be helpful
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param counts \code{logical}. Display counts of GC content rather than frequency
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param GCtheory \code{logical} default is \code{FALSE} to give the true GC content%, set to \code{TRUE} to normalize
#' values of GC_Content by the theoretical values using \code{\link{gcTheoretical}}. \code{species} must be specified.
#' @param GCtheoryType \code{"character"} Select type of data to normalize GC content% agianst accepts either "Genome" or
#' "Transcriptome" Default is "Genome"
#' @param GCobject an object of class GCTheoretical.
#'  Defaults to the gcTheoretical object supplied witht= the package
#' @param species \code{character} if \code{gcTheory} is \code{TRUE} its must be accompanied by a species
#' Species currently supported can be obtained using \code{mData(gcTheoretical)}
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
#' @return A plotly object
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
#'
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom viridisLite inferno
#' @importFrom magrittr %>%
#'
#' @export
plotGcHeatmap <- function(x, usePlotly = FALSE, counts = FALSE,
                          clusterNames = FALSE, pwfCols, labels,
                          GCtheory = FALSE, GCtheoryType = "Genome",
                          species = "Hsapiens", GCobject,
                          dendrogram = FALSE, ...){

  df <- tryCatch(Per_sequence_GC_content(x))

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Tidy up the GC content variables
  if(GCtheory & missing(GCobject)){
    GCobject <- ngsReports::gcTheoretical
  }

  if(GCtheory & GCtheoryType == "Genome"){
    spp <- genomes(GCobject)
    if(!species %in% spp$Name){
      stop(cat("Currently only supports genomes for", spp$Name, sep = ", "))
    }
  }

  if(GCtheory & GCtheoryType == "Transcriptome"){
    spp <- transcriptomes(GCobject)
    if(!species %in% spp$Name){
      stop(cat("Currently only supports transcriptomes for", spp$Name, sep = ", "))
    }
  }

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


  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- lapply(split(df, df$Filename), function(x){
      x$Value <-  x$Count / sum(x$Count)
      x
    }) %>%
      dplyr::bind_rows()

    if(GCtheory){
      df <- lapply(split(df, df$Filename), function(x){
        gcTheoryDF <- getGC(GCobject, name = species, type = GCtheoryType)
        x$Value <- x$Value - unlist(gcTheoryDF[species])
        x
      }) %>%
        dplyr::bind_rows()
    }
    fillLab <- "Frequency"
  }
  else{
    df$Value <- df$Count
    fillLab <- "Count"
  }

  if(clusterNames){
    # Grab the main columns & cast from long to wide
    mat <- reshape2::acast(df[c("Filename", "GC_Content", "Value")], Filename ~ GC_Content, value.var = "Value")
    mat[is.na(mat)] <- 0
    clus <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    mat <- mat[row.ord,]
    # Reform the data frame from the matrix
    df <- reshape2::melt(mat, varnames = c("Filename", "GC_Content"), value.name = "Value")
    df$Filename <- as.character(df$Filename)
    df$Value <- as.numeric(df$Value)
    df$GC_Content <- as.integer(df$GC_Content)
  }

  df$Filename <- labels[df$Filename]
  df$Filename <- factor(df$Filename, levels = unique(df$Filename))

  # Draw the heatmap
  GCheatmap <- ggplot(df, aes_string(x = "GC_Content", y = "Filename", fill = "Value")) +
    geom_tile() +
    scale_x_continuous(expand = c(0, 0)) +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank())

  if(GCtheory){
    GCheatmap <- GCheatmap +
      labs(fill = "Difference from\nTheoretical GC") +
      scale_fill_gradient2(low = inferno(1, begin = 0.4),
                           high = inferno(1, begin = 0.9),
                           midpoint = 0, mid = inferno(1, begin = 0))
  }
  else{
    GCheatmap <- GCheatmap +
      labs(fill = fillLab) +
      scale_fill_gradientn(colours = inferno(50))
  }
  if (!is.null(userTheme)) GCheatmap <- GCheatmap + userTheme

  if(usePlotly){

    GCheatmap <- GCheatmap +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

    t <- getSummary(x)
    t <- t[t$Category == "Per sequence GC content",]
    t$Filename <- factor(labels[t$Filename], levels = unique(df$Filename))
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
    t$x <- 1
    t$key <- labels[as.character(t$Filename)]

    if (missing(pwfCols)) pwfCols <- ngsReports::pwf
    col <- getColours(pwfCols)

    sideBar <- ggplot(t, aes_string(x = "x", y = "Filename", key = "key")) +
      geom_tile(aes_string(fill = "Status")) +
      scale_fill_manual(values = col) +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position="none",
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())
    sideBar <- suppressMessages(plotly::ggplotly(sideBar, tooltip = c("Status", "Filename")))

    #plot dendrogram
    if(dendrogram){
      ggdend <- function(df) {
        ggplot() +
          geom_segment(data = df, aes_string("x","y", xend = "xend", yend = "yend")) +
          ggdendro::theme_dendro()
      }

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 1)) +
        scale_x_continuous(expand = c(0,2))


      GCheatmap <- plotly::subplot(dendro, sideBar, GCheatmap, widths = c(0.1,0.1,0.8),
                                   margin = 0, shareY = TRUE) %>%
        plotly::layout(xaxis3 = list(title = "GC Content (%)"))
    }
    else{
      GCheatmap <- suppressWarnings(
        suppressMessages(
          plotly::subplot(plotly::plotly_empty(), sideBar, GCheatmap,
                          widths = c(0.1,0.1,0.8), margin = 0, shareY = TRUE) %>%
            plotly::layout(xaxis3 = list(title = "GC Content (%)"),
                           annotations = list(text = "Filename", showarrow = FALSE,
                                              textangle = -90))
        )
      )
    }

  }

  GCheatmap
}
