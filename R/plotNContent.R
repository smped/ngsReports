#' @title Draw an N Content Plot
#'
#' @description Draw an N Content Plot across one or more FastQC reports
#'
#' @details
#' This extracts the N_Content from the supplied object and generates a ggplot2
#' object, with a set of minimal defaults.
#' The output of this function can be further modified using the standard
#' ggplot2 methods.
#'
#' When `x` is a single FastqcData object line plots will always be drawn
#' for all Ns.
#' Otherwise, users can select line plots or heatmaps.
#'
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical`. Output as ggplot2 (default) or plotly
#' object.
#' @param warn,fail The default values for warn and fail are 5 and 10
#' respectively (i.e. percentages)
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default
#' @param lineCol Defaults to red
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param heat_w Relative width of any heatmap plot components
#' @param ... Used to pass additional attributes to theme() and between methods
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
#' # The default plot
#' plotNContent(fdl[[1]])
#'
#' @docType methods
#'
#' @import ggplot2
#' @importFrom dplyr vars
#' @import tibble
#' @importFrom tidyr unnest
#' @importFrom stringr str_split
#'
#' @name plotNContent
#' @rdname plotNContent-methods
#' @export
setGeneric("plotNContent", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
  standardGeneric("plotNContent")
}
)
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "ANY", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
  .errNotImp(x)
}
)
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "character", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
  x <- FastqcDataList(x)
  if (length(x) == 1) x <- x[[1]]
  plotNContent(x, usePlotly, labels, pwfCols, warn, fail, ...)
}
)
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...,
    lineCol = "red"){

  ## Get the NContent
  df <- getModule(x, "Per_base_N_content")
  colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

  ## Handle empty/missing modules
  msg <- c()
  if (!length(df)) msg <- "No N Content Module Detected"
  if (sum(df[["Percentage"]]) == 0) msg <- "No N Content in Sequences"
  if (!is.null(msg)) {
    nPlot <- ngsReports:::.emptyPlot(msg)
    if (usePlotly) nPlot <- ggplotly(nPlot, tooltip = "")
    return(nPlot)
  }

  ## Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(.isValidPwf(pwfCols))
  pwfCols <- setAlpha(pwfCols, 0.2)

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]
  df$Base <- factor(df$Base, levels = unique(df$Base))
  df$xValue <- as.integer(df$Base)

  ## Setup the BG colours
  rects <- tibble(
    xmin = 0,
    xmax = max(df$xValue),
    ymin = c(0, warn, fail),
    ymax = c(warn, fail, 100),
    Status = c("PASS", "WARN", "FAIL")
  )

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  yLab <- "N Content (%)"
  x <- "xValue"
  nPlot <- ggplot(df) +
    geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
    ) +
    geom_line(aes(!!sym(x), Percentage), colour = lineCol) +
    geom_point(
      aes(!!sym(x), Percentage, group = Base),
      size = 0, colour = rgb(0, 0, 0, 0)
    ) +
    scale_fill_manual(values = getColours(pwfCols)) +
    scale_x_continuous(
      breaks = unique(df$xValue), labels = levels(df$Base), expand = c(0,0)
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    facet_wrap(~Filename) +
    labs(x = "Position in Read (bp)", y = yLab) +
    guides(fill = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ## Add the basic customisations
  if (!is.null(userTheme)) nPlot <- nPlot + userTheme

  if (usePlotly) {
    nPlot <- nPlot +
      xlab("") +
      theme(legend.position = "none")
    nPlot <- suppressMessages(plotly::ggplotly(nPlot))
    nPlot <- suppressMessages(
      suppressWarnings(
        plotly::subplot(
          plotly::plotly_empty(),
          nPlot,
          widths = c(0.14,0.86)
        )
      )
    )
    nPlot <- plotly::layout(nPlot, yaxis2 = list(title = yLab))


    ## Set the hoverinfo for bg rectangles to the vertices only,
    ## This will effectively hide them
    nPlot$x$data <- lapply(nPlot$x$data, .hidePWFRects)
    ## Hide the xValue parameter to make it look nicer
    nPlot$x$data[[6]]$text <- gsub(
      "(.+)(xValue.+)(Percentage.+)",
      "\\1\\3",
      nPlot$x$data[[6]]$text
    )
  }
  nPlot
}
)
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20,
    cluster = FALSE, dendrogram = FALSE, heat_w = 8, ...){

  ## Get the NContent
  df <- getModule(x, "Per_base_N_content")
  colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

  ## Handle empty/missing modules
  msg <- c()
  if (!length(df)) msg <- "No N Content Module Detected"
  if (sum(df[["Percentage"]]) == 0) msg <- "No N Content in Sequences"
  if (!is.null(msg)) {
    nPlot <- .emptyPlot(msg)
    if (usePlotly) nPlot <- ggplotly(nPlot, tooltip = "")
    return(nPlot)
  }

  ## Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(.isValidPwf(pwfCols))

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]

  ## fill bins up to the max sequence length
  df$Base <- lapply(
    df$Base,
    function(x){
      rng <- as.integer(str_split(x, pattern = "-")[[1]])
      seq(min(rng), max(rng), by = 1L)
    }
  )
  df <- unnest(df, Base)

  ## Now define the order for a dendrogram if required
  key <- names(labels)
  cols <- c("Filename", "Base", "Percentage")
  clusterDend <- .makeDendro(df[cols], "Filename", "Base", "Percentage")
  dx <- ggdendro::dendro_data(clusterDend)
  if (dendrogram | cluster) key <- labels(clusterDend)
  if (!dendrogram) dx$segments <- dx$segments[0,]
  df$Filename <- factor(labels[df$Filename], levels = labels[key])

  cols <- .makePwfGradient(
    df$Percentage, pwfCols,
    breaks = c(0, warn, fail, 101), passLow = TRUE,
    na.value = "white"
  )

  xLab <- "Position in Read (bp)"
  hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
  nPlot <- ggplot(df, aes(Base, Filename, fill = Percentage, label = Base)) +
    geom_tile() +
    do.call("scale_fill_gradientn", cols) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), position = "right") +
    labs(x = xLab, y = NULL, fill = "%N") +
    ggtitle("Per Base N Content") +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"),
      plot.title = element_text(hjust = hj),
      axis.title.x = element_text(hjust = hj)
    )

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
  if (!is.null(userTheme)) nPlot <- nPlot + userTheme

  ## Reset the status using current values
  status <- dplyr::summarise(
    dplyr::group_by(df, Filename),
    Percentage = max(Percentage, na.rm = TRUE)
  )
  status$Status <- cut(
    status$Percentage,
    breaks = c(0, warn, fail, 101),
    include.lowest = TRUE,
    labels = c("PASS", "WARN", "FAIL")
  )
  status$Filename <- factor(labels[status$Filename], levels = labels[key])

  .prepHeatmap(nPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

}
)
