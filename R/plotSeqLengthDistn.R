#' @title Plot the Sequence Length Distribution
#'
#' @description Plot the Sequence Length Distribution across one or more FASTQC
#' reports
#'
#' @details
#' This extracts the Sequence Length Distribution from the supplied object and
#' generates a ggplot2 object, with a set of minimal defaults.
#' The output of this function can be further modified using the standard
#' ggplot2 methods.
#'
#' A cdf plot can also be generated to provide guidance for minimum
#' read length in some NGS workflows, by setting `plotType = "cdf"`.
#' If all libraries have reads of identical lengths, these plots may be less
#' informative.
#'
#' An alternative interactive plot is available by setting the argument
#' `usePlotly = TRUE`.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical`. Output as ggplot2 or plotly object.
#' @param plotType `character`. Can only take the values
#' `plotType = "heatmap"` `plotType = "line"` or
#' `plotType = "cdf"`
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param counts `logical` Should distributions be shown as counts or
#' frequencies (percentages)
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` and
#' `usePlotly` are `FALSE`. If both `cluster` and
#' `dendrogram` are specified as `TRUE` then the dendrogram
#' will be displayed.
#' @param heat_w Relative width of any heatmap plot components
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail values in plot
#' @param ... Used to pass additional attributes to theme()
#' @param expand.x Output from `expansion()` or numeric vector of
#' length 4. Passed to `scale_x_discrete`
#' @param heatCol The colour scheme for the heatmap
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Plot as a frequency plot using lines
#' plotSeqLengthDistn(fdl)
#'
#' # Or plot the cdf
#' plotSeqLengthDistn(fdl, plotType = "cdf")
#'
#' @docType methods
#'
#' @importFrom dplyr vars mutate across
#' @importFrom plotly ggplotly layout subplot
#' @importFrom tidyr pivot_wider complete nesting unnest
#' @importFrom tidyselect everything any_of
#' @importFrom grDevices hcl.colors
#' @import ggplot2
#'
#' @name plotSeqLengthDistn
#' @rdname plotSeqLengthDistn-methods
#' @export
setGeneric("plotSeqLengthDistn", function(
    x, usePlotly = FALSE, labels, ...){
  standardGeneric("plotSeqLengthDistn")
}
)
#' @rdname plotSeqLengthDistn-methods
#' @export
setMethod("plotSeqLengthDistn", signature = "ANY", function(
    x, usePlotly = FALSE, labels, ...){
  .errNotImp(x)
}
)
#' @rdname plotSeqLengthDistn-methods
#' @export
setMethod("plotSeqLengthDistn", signature = "character", function(
    x, usePlotly = FALSE, labels, ...){
  x <- FastqcDataList(x)
  if (length(x) == 1) x <- x[[1]]
  plotSeqLengthDistn(x, usePlotly, labels, ...)
}
)
#' @rdname plotSeqLengthDistn-methods
#' @export
setMethod("plotSeqLengthDistn", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, plotType = c("line", "cdf"), ...,
    expand.x = expansion(0, 0.2)){

  df <- getModule(x, "Sequence_Length_Distribution")
  plotType <- match.arg(plotType)

  if (!length(df)) {
    lenPlot <- .emptyPlot("No Sequence Length Module Detected")
    if (usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
    return(lenPlot)
  }

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]

  ## Add zero counts for lengths either side of the included range
  ## This is only required if a single value exists
  if (nrow(df) == 1) {
    df <- dplyr::bind_rows(
      df,
      dplyr::mutate(df, Lower = Lower - 1, Count = 0),
      dplyr::mutate(df, Lower = Lower + 1, Count = 0)
    )
  }

  df$Lower <- as.integer(df$Lower)
  df <- dplyr::arrange_at(df, vars("Lower"))
  df <- df[c("Filename", "Length", "Lower", "Count")]
  df$Cumulative <- cumsum(df$Count)
  df$Length <- factor(df$Lower, levels = unique(df$Lower))

  ## Sort out some plotting parameters
  stopifnot(is.numeric(expand.x), length(expand.x) == 4)
  xLab <- "Sequence Length (bp)"
  yLab <- c(cdf = "Cumulative Count", line = "Count")[plotType]
  plotY <- c(cdf = "Cumulative", line = "Count")[plotType]

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  lenPlot <- ggplot(
    df,
    aes_string("Length", plotY, colour = "Filename", group = "Filename")
  ) +
    geom_line() +
    facet_wrap(~Filename) +
    labs(x = xLab, y = yLab) +
    scale_x_discrete(expand = expand.x) +
    scale_y_continuous(labels = scales::comma) +
    theme_bw() +
    theme(legend.position = "none")

  if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

  if (usePlotly) {
    lenPlot <-
      suppressMessages(plotly::ggplotly(lenPlot, tooltip = c("x", "y")))

    lenPlot <- suppressMessages(
      suppressWarnings(
        plotly::subplot(
          plotly::plotly_empty(),
          lenPlot,
          widths = c(0.14,0.86))
      ))
    lenPlot <- plotly::layout(
      lenPlot,
      xaxis2 = list(title = xLab),
      yaxis2 = list(title = yLab)
    )

  }

  lenPlot

}
)
#' @rdname plotSeqLengthDistn-methods
#' @export
setMethod(
  "plotSeqLengthDistn", signature = "FastqcDataList",
  function(
    x, usePlotly = FALSE, labels, counts = FALSE,
    plotType = c("heatmap", "line", "cdf"), cluster = FALSE,
    dendrogram = FALSE, heat_w = 8, pwfCols, ...,
    heatCol = hcl.colors(50, "inferno")
  ){

    mod <- "Sequence_Length_Distribution"
    df <- getModule(x, mod)
    if (!length(df)) {
      lenPlot <- .emptyPlot("No Sequence Length Module Detected")
      if (usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
      return(lenPlot)
    }

    ## Check for valid plotType
    plotType <- match.arg(plotType)

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(x, labels, ...)
    labels <- labels[names(labels) %in% df$Filename]

    ## Lengths will probably be binned so define the bins then expand
    ## the range at the lower and upper limits to add zero.
    ## This will enable replication of the default FastQC plot
    ## In reality, this will only be required when there are <2 bins
    df <- complete(
      df, Filename, nesting(Length, Lower, Upper), fill = list(Count = 0)
    )
    df <- dplyr::arrange(df, Lower)
    df$Length <- forcats::fct_inorder(df$Length)
    df <- dplyr::arrange(df, Filename, Length)

    ## Get the cdf count
    df <- dplyr::group_by(df, Filename)
    df <- mutate(df, Cumulative = cumsum(Count))
    if (!counts) {
      df <- mutate(
        df, Cumulative = Cumulative / max(Cumulative), Freq = Count / sum(Count)
      )
    }
    df <- dplyr::ungroup(df)
    ## Round the values for better plotting
    df <- mutate(df, across(any_of(c("Cumulative", "Freq")), round, digits = 4))

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    ## Check axis expansion
    if (missing(pwfCols)) pwfCols <- pwf

    if (plotType %in% c("line", "cdf")) {

      ## Decide whether to plot the Count or cdf sum
      ## and set all labels
      plotY <- dplyr::case_when(
        plotType == "cdf" ~ "Cumulative",
        plotType == "line" & counts ~ "Count",
        plotType == "line" & !counts ~ "Freq"
      )
      yLab <- dplyr::case_when(
        plotType == "cdf" & counts ~ "Cumulative Count",
        plotType == "cdf" & !counts ~ "Cumulative (%)",
        plotType == "line" & counts ~ "Count",
        plotType == "line" & !counts ~ "Percent (%)"
      )
      yLabelFun <- ifelse(counts, scales::comma, scales::percent)

      df$Filename <- labels[df$Filename]
      lenPlot <- ggplot(
        df,
        aes_string(
          x = "Length", y = plotY, colour = "Filename", group = "Filename"
        )
      ) +
        geom_line() +
        labs(y = yLab) +
        scale_y_continuous(labels = yLabelFun) +
        theme_bw()

      if (usePlotly) {
        ttip <- c("x", "y", "colour")
        lenPlot <- suppressMessages(
          suppressWarnings(
            plotly::ggplotly(
              lenPlot + theme(legend.position = "none"), tooltip = ttip
            )
          )
        )
      }
    }

    if (plotType == "heatmap") {

      ## Now define the order for a dendrogram if required
      ## This only applies to a heatmap
      key <- names(labels)
      cols <- c("Filename", "Length", "Freq")
      clusterDend <- .makeDendro(df[cols], "Filename","Length", "Freq")
      dx <- ggdendro::dendro_data(clusterDend)
      if (dendrogram | cluster) key <- labels(clusterDend)
      df$Filename <- factor(labels[df$Filename], levels = labels[key])
      if (!dendrogram) dx$segments <- dx$segments[0,]
      df$y <- as.integer(df$Filename)
      df$Percent <- scales::percent(df$Freq, accuracy = 0.1)
      df$Total <- scales::comma(df$Count)

      ## Make the basic heatmap. The first aes sets the labels for plotly
      aes <- aes(
        `%` = Percent, Count = Total, Length = Length, Filenane = Filename
      )
      hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
      lenPlot <- ggplot(df, aes) +
        geom_rect(
          aes_string(
            xmin = "Lower", xmax = "Upper + 1",
            ymin = "y - 0.5", ymax = "y + 0.5", fill = "Freq"
          ),
          colour = NA
        ) +
        labs(x = "Sequence Length", fill = "Percent") +
        ggtitle(gsub("_", " ", mod)) +
        scale_fill_gradientn(
          colours = heatCol, labels = scales::percent, limits = c(0, 1)
        ) +
        scale_colour_gradientn(
          colours = heatCol, labels = scales::percent, limits = c(0, 1)
        ) +
        guides(colour = "none") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(
          expand = c(0, 0), position = "right",
          breaks = seq_along(levels(df$Filename)), labels = levels(df$Filename)
        ) +
        theme_bw() +
        theme(
          plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"),
          plot.title = element_text(hjust = hj),
          axis.title.x = element_text(hjust = hj)
        )
      if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

      ## Get the PWF status
      status <- subset(getSummary(x), Category == gsub("_", " ", mod))
      status$Filename <- factor(labels[status$Filename], levels = labels[key])

      lenPlot <- .prepHeatmap(
        lenPlot, status, dx$segments, usePlotly, heat_w, pwfCols
      )

    }
    lenPlot
  }
)
