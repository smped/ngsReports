#' @title Plot the Per Sequence Quality Scores
#'
#' @description Plot the Per Sequence Quality Scores for a set of FASTQC reports
#'
#' @details Plots the distribution of average sequence quality scores across the
#' set of files. Values can be plotted either as counts (`counts = TRUE`)
#' or as frequencies (`counts = FALSE`).
#'
#' Any faceting or scale adjustment can be performed after generation of the
#' initial plot, using the standard methods of ggplot2 as desired.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or path
#' @param counts `logical`. Plot the counts from each file if
#' `counts = TRUE`, otherwise the frequencies will be plotted
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType `character`. Can only take the values
#' `plotType = "heatmap"` or `plotType = "line"`
#' @param warn,fail The default values for warn and fail are 5 and 10
#' respectively (i.e. percentages)
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param alpha set alpha for line graph bounds
#' @param heatCols Colour palette for the heatmap
#' @param heat_w Relative width of any heatmap plot components
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
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
#' plotSeqQuals(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fqName(fdl))
#' plotSeqQuals(fdl[r1])
#'
#' @docType methods
#'
#' @importFrom stats hclust dist
#' @importFrom scales percent comma percent_format
#' @importFrom tidyr pivot_wider complete nesting
#' @importFrom grDevices hcl.colors
#' @importFrom rlang "!!" sym
#' @import ggplot2
#'
#' @name plotSeqQuals
#' @rdname plotSeqQuals-methods
#' @export
setGeneric("plotSeqQuals", function(
    x, usePlotly = FALSE, labels, pwfCols, counts = FALSE, alpha = 0.1,
    warn = 30, fail = 20, ...){
  standardGeneric("plotSeqQuals")
}
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod("plotSeqQuals", signature = "ANY", function(
    x, usePlotly = FALSE, labels, pwfCols, counts = FALSE, alpha = 0.1,
    warn = 30, fail = 20, ...){
  .errNotImp(x)
}
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod("plotSeqQuals", signature = "character", function(
    x, usePlotly = FALSE, labels, pwfCols, counts = FALSE, alpha = 0.1,
    warn = 30, fail = 20, ...){
  x <- FastqcDataList(x)
  if (length(x) == 1) x <- x[[1]]
  plotSeqQuals(
    x, usePlotly, labels, pwfCols, counts, alpha, warn, fail, ...)
}
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod("plotSeqQuals", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, pwfCols, counts = FALSE, alpha = 0.1,
    warn = 30, fail = 20, ...){

  df <- getModule(x, "Per_sequence_quality_scores")

  if (!length(df)) {
    qualPlot <- .emptyPlot("No Sequence Quality Moudule Detected")
    if(usePlotly) qualPlot <- ggplotly(qualPlot, tooltip = "")
    return(qualPlot)
  }

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]

  ## Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(.isValidPwf(pwfCols))
  pwfCols <- setAlpha(pwfCols, alpha)
  stopifnot(warn > fail)

  ## Find the minimum quality value
  minQ <- min(df$Quality)

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  ## make Ranges for rectangles and set alpha
  rects <- tibble(
    ymin = 0,
    ymax = max(df$Count),
    xmin = c(0, fail, warn),
    xmax = c(fail, warn, 41),
    Status = c("FAIL", "WARN", "PASS")
  )

  xLab <- "Mean Sequence Quality Per Read (Phred Score)"
  yLab <- "Number of Sequences"

  if (!counts){

    ## Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename)
    df <- dplyr::mutate(df, Frequency = Count / sum(Count))
    df <- dplyr::ungroup(df)
    df$Frequency <- round(df$Frequency, 3)
    df$Percent <- scales::percent(df$Frequency)
    rects$ymax <- max(df$Frequency)

    qualPlot <- ggplot(df) +
      geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
      ) +
      geom_line(aes(Quality, Frequency, colour = Filename))
    yLabelFun <- scales::percent_format(accuracy = 1)
    yLab <- "Frequency"

  }
  else{
    ## Initialise the plot using counts
    qualPlot <- ggplot(df) +
      geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
      ) +
      geom_line(aes(x = Quality, y = Count, colour = Filename))
    yLabelFun <- scales::comma

  }

  qualPlot <- qualPlot +
    scale_fill_manual(values = getColours(pwfCols))  +
    scale_y_continuous(
      limits = c(0, rects$ymax[1]), expand = c(0, 0), labels = yLabelFun
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_colour_discrete()  +
    facet_wrap(~Filename) +
    labs(x = xLab, y = yLab) +
    theme_bw() +
    theme(legend.position = "none")

  if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

  if (usePlotly) {

    ## Render as a plotly object
    hv <-  c("x", "y", "colour")
    qualPlot <- suppressMessages(
      suppressWarnings(plotly::ggplotly(qualPlot, hoverinfo = hv))
    )

    qualPlot <- suppressMessages(
      suppressWarnings(
        plotly::subplot(
          plotly::plotly_empty(), qualPlot, widths = c(0.14,0.86)
        )
      )
    )
    qualPlot <- plotly::layout(
      qualPlot, xaxis2 = list(title = xLab), yaxis2 = list(title = yLab)
    )

    ## Set the hoverinfo for bg rectangles to the vertices only,
    ## This will effectively hide them
    qualPlot$x$data <- lapply(qualPlot$x$data, function(x){
      ## If there is a name component & it contains
      ## PASS/WARN/FAIL set the hoverinfo to none
      if ("name" %in% names(x)) {
        if (grepl("(PASS|WARN|FAIL)", x$name)) {
          x$hoverinfo <- "none"
        }
      }
      x
    })
  }

  ## Draw the plot
  qualPlot

}
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod("plotSeqQuals", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, pwfCols, counts = FALSE, alpha = 0.1,
    warn = 30, fail = 20, plotType = c("heatmap", "line"), dendrogram = FALSE,
    cluster = FALSE, heatCols = hcl.colors(100, "inferno"), heat_w = 8, ...){

  ## Read in data
  mod <- "Per_sequence_quality_scores"
  df <- getModule(x, "Per_sequence_quality_scores")

  if (!length(df)) {
    qualPlot <- .emptyPlot("No Sequence Quality Moudule Detected")
    if (usePlotly) qualPlot <- ggplotly(qualPlot, tooltip = "")
    return(qualPlot)
  }

  ## As the Quality Scores may have different ranges across different
  ## samples, fill missing values with zero
  df <- complete(df, Filename, nesting(Quality), fill = list(Count = 0))
  df$Quality <- as.integer(df$Quality)

  ## Check for valid plotType
  plotType <- match.arg(plotType)
  xLab <- "Mean Sequence Quality Per Read (Phred Score)"
  plotVal <- ifelse(counts, "Count", "Frequency")
  stopifnot(warn > fail)

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]

  ## Sort out the colours
  if (base::missing(pwfCols)) pwfCols <- pwf

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  if (plotType == "heatmap") {

    ## Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename)
    df <- dplyr::mutate(df, Frequency = Count / sum(Count))
    df <- dplyr::ungroup(df)
    df$Frequency <- round(df$Frequency, 3)
    df$Total <- scales::comma(df$Count)

    ## Now define the order for a dendrogram if required
    ## This only applies to a heatmap
    key <- names(labels)
    cols <- c("Filename", "Quality", "Frequency")
    clusterDend <- .makeDendro(df[cols], "Filename","Quality", "Frequency")
    dx <- ggdendro::dendro_data(clusterDend)
    if (dendrogram | cluster) key <- labels(clusterDend)
    df$Filename <- factor(labels[df$Filename], levels = labels[key])
    if (!dendrogram) dx$segments <- dx$segments[0,]

    hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
    qualPlot <- ggplot(df, aes(Quality, Filename, C = Total)) +
      geom_tile(aes(fill = Frequency)) +
      ggtitle("Per Sequence Quality Score") +
      labs(x = xLab, y = c()) +
      scale_fill_gradientn(colours = heatCols, limits = c(0, 1)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0), position = "right") +
      theme(
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0), "points"),
        plot.title = element_text(hjust = hj),
        axis.title.x = element_text(hjust = hj)
      )

    status <- getSummary(x)
    status <- subset(status, Category == gsub("_", " ", mod))
    status$Filename <- factor(labels[status$Filename], levels = labels[key])

    qualPlot <-
      .prepHeatmap(qualPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

  }

  if (plotType == "line") {

    ## make Ranges for rectangles and set alpha
    pwfCols <- setAlpha(pwfCols, alpha)
    rects <- tibble(
      ymin = 0,
      ymax = max(df$Count),
      xmin = c(0, fail, warn),
      xmax = c(fail, warn, 41),
      Status = c("FAIL", "WARN", "PASS")
    )
    ## No clustering required so just use the labels
    df$Filename <- labels[df$Filename]

    yLab <- ifelse(counts, "Number of Sequences", "Frequency")
    yLabelFun <-
      ifelse(counts, scales::comma, scales::percent_format(accuracy = 1))
    plotVal <- ifelse(counts, "Count", "Frequency")

    if (!counts) {

      ## Summarise to frequencies & initialise the plot
      df <- dplyr::group_by(df, Filename)
      df <- dplyr::mutate(df, Frequency = Count / sum(Count))
      df <- dplyr::ungroup(df)
      df$Frequency <- round(df$Frequency, 4)
      rects$ymax <- max(df$Frequency)

    }

    qualPlot <- ggplot(df) +
      geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
      ) +
      geom_line(aes(x = Quality, y = !!sym(plotVal), colour = Filename)) +
      scale_fill_manual(values = getColours(pwfCols))  +
      scale_y_continuous(
        limits = c(0, rects$ymax[1]), expand = c(0, 0), labels = yLabelFun
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_colour_discrete() +
      labs(x = xLab, y = yLab) +
      guides(fill = "none") +
      theme_bw()

    if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

    if (usePlotly) {

      hv <- c("x", "y", "colour")
      qualPlot <- qualPlot + theme(legend.position = "none")
      qualPlot <- suppressMessages(
        suppressWarnings(
          plotly::ggplotly(qualPlot, hoverinfo = hv)
        )
      )

      ## Turn off the hoverinfo for the bg rectangles
      ## This will effectively hide them
      qualPlot$x$data <- lapply(qualPlot$x$data, function(x){
        ## If there is a name component & it contains
        ## PASS/WARN/FAIL set the hoverinfo to none
        if ("name" %in% names(x)) {
          if (grepl("(PASS|WARN|FAIL)", x$name)) {
            x$hoverinfo <- "none"
          }
        }
        x
      })
    }}

  qualPlot
}

)
