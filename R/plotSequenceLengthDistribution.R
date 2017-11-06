#' @title Draw an Adapter Content Plot
#'
#' @description Draw an Adapter Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the Adapter_Content from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#'
#' Preset axis limits can also be overwritten easily by adding a call to \code{scale_y_continuous}
#' after the call to \code{plotAdapterContent}.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical}. Output as ggplot2 or plotly object.
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param counts \code{logical} Should distributions be shown as counts or frequencies (percentages)
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} and \code{usePlotly} are \code{FALSE}.
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param ... Used to pass additional attributes to theme()
#' @param expand.x Passed to \code{scale_x_discrete}
#' @param lineWidth,lineCol Passed to geom_hline and geom_vline to determine
#' width and colour of gridlines
#' @param heatCol The colour scheme for the heatmap
#'
#' @return A standard ggplot2 object, or an interactive plotly object
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
#' # Plot as a frequency plot using lines
#' plotSequenceLengthDistribution(fdl)
#'
#' @importFrom dplyr vars
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom plotly ggplotly
#' @importFrom plotly layout
#' @importFrom plotly subplot
#' @importFrom viridisLite inferno
#'
#' @export
plotSequenceLengthDistribution <- function(x, usePlotly = FALSE, labels, counts = FALSE,
                                           plotType = "heatmap", clusterNames = FALSE,
                                           dendrogram = FALSE, ...,
                                           expand.x = c(0,0.2), lineCol = "grey20",
                                           lineWidth = 0.2, heatCol = inferno(50)){

  df <- tryCatch(Sequence_Length_Distribution(x))

  # Check for valid plotType
  stopifnot(plotType %in% c("line", "heatmap"))

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Add zero counts for lengths either side of the included range
  df <- dplyr::bind_rows(
    lapply(split(df, f = df$Filename),
           function(x){
             # Fix global environment error
             Lower <- NULL
             dplyr::bind_rows(x,
                              dplyr::data_frame(Filename = x$Filename[1],
                                                Lower = c(min(x$Lower) - 1, max(x$Upper) + 1),
                                                Upper = Lower,
                                                Length = as.character(Lower),
                                                Count = 0)
             )
           })
  )

  df$Lower <- as.integer(df$Lower)
  df <- reshape2::dcast(df, Filename ~ Lower, value.var = "Count")
  df[is.na(df)] <- 0


  if(clusterNames){
    xx <- dplyr::select(df, -Filename)
    clus <- as.dendrogram(hclust(dist(xx)))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
  }
  df$Filename <- labels[df$Filename]
  df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Lower", value.name = "Count")



  df$Filename <- factor(df$Filename, levels = unique(df$Filename))


  df$Count <- as.numeric(df$Count)
  # Convert the counts to frequencies
  df <- dplyr::bind_rows(
    lapply(split(df, f = df$Filename),
           function(x){
             x$Freq <- 100*x$Count / sum(x$Count)
             x
           }))

  # Arrange in position
  df <- dplyr::arrange_at(df, vars("Filename", "Lower"))
  df$Lower <- factor(df$Lower, levels = unique(df$Lower))

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  if (plotType == "line"){
    if (counts){
      lenPlot <- ggplot(df, aes_string("Lower", "Count", colour = "Filename", group = "Filename")) +
        geom_line() +
        labs(y = "Count") +
        theme_bw()
    }
    else{
      lenPlot <- ggplot(df, aes_string("Lower", "Freq", colour = "Filename", group = "Filename")) +
        geom_line() +
        scale_y_continuous() +
        labs(y = "Percent (%)") +
        theme_bw()
    }
  }

  if (plotType == "heatmap"){

    if (counts){
      lenPlot <- ggplot(df, aes_string("Lower","Filename", fill = "Count")) +
        geom_tile(colour = lineCol) +
        scale_fill_gradientn(colours = heatCol) +
        scale_y_discrete(labels = labels, expand = c(0, 0))
    }
    else{
      lenPlot <- ggplot(df, aes_string(x = "Lower",y ="Filename", fill = "Freq")) +
        geom_tile(colour = lineCol) +
        labs(fill = "Percent (%)") +
        scale_fill_gradientn(colours = heatCol) +
        scale_y_discrete(labels = labels, expand = c(0, 0))
    }
  }

  lenPlot <- lenPlot +
    scale_x_discrete(expand = expand.x) +
    labs(x = "Sequence Length") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

  if(usePlotly){
    if(plotType == "line") {
      lenPlot <- lenPlot +
        theme(legend.position = "none")
      lenPlot <- ggplotly(lenPlot, tooltip = c("x", "y", "colour"))
    }else{

      lenPlot <- lenPlot  +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank())

      pwfCols <- ngsReports::pwf

      lenPlot <- lenPlot +
        theme(panel.background = element_blank())

      t <- getSummary(x)
      t <- t[t$Category == "Sequence Length Distribution",]
      t$Filename <- labels[t$Filename]
      t$Filename <- factor(t$Filename, levels = levels(df$Filename))
      t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")

      # Set the key for extracting individual files on click
      key <- names(labels[match(levels(df$Filename), labels)])

     sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)


      #plot dendrogram
      if(dendrogram && clusterNames){

        dx <- ggdendro::dendro_data(clus)
        dendro <- ggdend(dx$segments) +
          coord_flip() +
          scale_y_reverse(expand = c(0, 0)) +
          scale_x_continuous(expand = c(0,0.5))
      }
      else{
        dendro <- plotly::plotly_empty()
      }

      lenPlot <- suppressMessages(
        plotly::subplot(dendro, sideBar, lenPlot, widths = c(0.1, 0.08, 0.82),
                        margin = 0.001, shareY = TRUE)
      )
    }

  }

  lenPlot

}
