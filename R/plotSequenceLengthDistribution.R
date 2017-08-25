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
#' @param ... Used to pass additional attributes to theme()
#' @param expand.x Passed to \code{scale_x_discrete}
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
#' @importFrom ggplot2 scale_colour_discrete
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#'
#' @export
plotSequenceLengthDistribution <- function(x, usePlotly = FALSE, labels, counts = FALSE,
                                           plotType = "heatmap",
                                           ...,
                                           expand.x = c(0,0.2)){

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
    lapply(split(df, f = df$Filename) ,
           function(x){
             dplyr::bind_rows(x,
                              data_frame(Filename = x$Filename,
                                         Lower = c(min(x$Lower) - 1, max(x$Upper) + 1),
                                         Upper = Lower,
                                         Length = as.character(Lower),
                                         Count = 0)
             )
           })
  )

  # Convert the counts to frequencies
  df <- dplyr::bind_rows(
    lapply(split(df, f = df$Filename),
           function(x){
             x$Freq <- 100*x$Count / sum(x$Count)
             x
           }))

  # Arrange in position
  df <- dplyr::arrange_at(df, vars("Filename", "Lower"))
  df$Length <- factor(df$Length, levels = unique(df$Length))

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  if (plotType == "line"){
    if (counts){
      lenPlot <- ggplot(df, aes_string("Length", "Count", colour = "Filename", group = "Filename")) +
        geom_line() +
        labs(y = "Count")
    }
    else{
      lenPlot <- ggplot(df, aes_string("Length", "Freq", colour = "Filename", group = "Filename")) +
        geom_line() +
        scale_y_continuous(limits = c(0, 100)) +
        labs(y = "Percent (%)")
    }
  }

  if (plotType == "heatmap"){

    df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

    if (counts){
      lenPlot <- ggplot(df, aes_string("Length","Filename", fill = "Count")) +
        geom_tile() +
        scale_fill_gradientn(colours = viridisLite::inferno(50)) +
        scale_y_discrete(labels = labels, expand = c(0, 0))
    }
    else{
      lenPlot <- ggplot(df, aes_string("Length","Filename", fill = "Freq")) +
        geom_tile() +
        labs(fill = "Percent (%)") +
        scale_fill_gradientn(colours = viridisLite::inferno(50)) +
        scale_y_discrete(labels = labels, expand = c(0, 0))
    }
  }

  lenPlot <- lenPlot +
    scale_x_discrete(expand = expand.x) +
    labs(x = "Sequence Length") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

  lenPlot

}
