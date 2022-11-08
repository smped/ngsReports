#' @title Plot the Base Qualities for each file
#'
#' @description Plot the Base Qualities for each file as separate plots
#'
#' @details When acting on a `FastqcDataList`, this defaults to a heatmap
#' using the mean Per_base_sequence_quality score. A set of plots which
#' replicate those obtained through a standard FastQC html report can be
#' obtained by setting `plotType = "boxplot"`, which uses `facet_wrap`
#' to provide the layout as a single ggplot object.
#'
#' When acting an a `FastqcData` object, this replicates the
#' `Per base sequence quality` plots from FastQC with no faceting.
#'
#' For large datasets, subsetting by R1 or R2 reads may be helpful.
#'
#' An interactive plot can be obtained by setting `usePlotly = TRUE`.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or character
#' vector of file paths
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param nc `numeric`. The number of columns to create in the plot layout.
#' Only used if drawing boxplots for multiple files in a FastqcDataList
#' @param warn,fail The default values for warn and fail are 30 and 20
#' respectively (i.e. percentages)
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType `character` Can be either `"boxplot"` or
#' `"heatmap"`
#' @param plotValue `character` Type of data to be presented. Can be
#' any of the columns returned by
#' `getModule(x, module = "Per_base_sequence_qual")`
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail values in plot
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param boxWidth set the width of boxes when using a boxplot
#' @param heat_w Relative width of any heatmap plot components
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A standard ggplot2 object or an interactive plotly object
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
#' # The default plot for multiple libraries is a heatmap
#' plotBaseQuals(fdl)
#'
#' # The default plot for a single library is the standard boxplot
#' plotBaseQuals(fdl[[1]])
#'
#' @docType methods
#'
#' @import ggplot2
#' @importFrom stats as.dendrogram order.dendrogram na.omit hclust dist
#' @import tibble
#' @importFrom tidyr unnest
#' @importFrom stringr str_split
#' @importFrom rlang "!!" sym
#'
#' @name plotBaseQuals
#' @rdname plotBaseQuals-methods
#' @export
setGeneric("plotBaseQuals", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 25, fail = 20,
    boxWidth = 0.8, ...){
  standardGeneric("plotBaseQuals")
}
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod("plotBaseQuals", signature = "ANY", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 25, fail = 20,
    boxWidth = 0.8, ...){
  .errNotImp(x)
}
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod("plotBaseQuals", signature = "character", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 25, fail = 20,
    boxWidth = 0.8, ...){
  x <- FastqcDataList(x)
  if (length(x) == 1) x <- x[[1]]
  plotBaseQuals(x, usePlotly, labels, pwfCols, warn, fail, boxWidth, ...)
}
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod("plotBaseQuals", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 25, fail = 20,
    boxWidth = 0.8, ...){

  ## Get the data
  df <- getModule(x, "Per_base_sequence_quality")

  ## Make a blank plot if no data is found
  if (!length(df)) {
    msg <- "No Per Base Sequence Quality Module Detected"
    qualPlot <- .emptyPlot(msg)
    if (usePlotly) qualPlot <- ggplotly(qualPlot, tooltip = "")
    return(qualPlot)
  }

  # Get the plot labels organised
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]

  stopifnot(is.numeric(boxWidth))
  df$Base <- factor(df$Base, levels = unique(df$Base))
  df$Position <- as.integer(df$Base)
  df$xmin <- df$Position - boxWidth/2
  df$xmax <- df$Position + boxWidth/2

  ## Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(.isValidPwf(pwfCols))
  pwfCols <- setAlpha(pwfCols, 0.2)

  ## Get any theme arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  ## Set the limits & rectangles
  ylim <- c(0, max(df$`90th_Percentile`) + 1)
  expand_x <- round(0.015*(max(df$Position) - min(df$Position)), 1)
  rects <- tibble(
    xmin = min(df$Position) - expand_x,
    xmax = max(df$Position) + expand_x,
    ymin = c(0, fail, warn),
    ymax = c(fail, warn, max(ylim)),
    Status = c("FAIL", "WARN", "PASS")
  )

  ## Get the Illumina encoding
  enc <- getModule(x, "Basic_Statistics")$Encoding[1]
  enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)
  ylab <- paste0("Quality Scores (", enc, " encoding)")

  ## Generate the basic plot
  lq <- "Lower_Quartile"
  uq <- "Upper_Quartile"
  qualPlot <- ggplot(df) +
    geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
    ) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = !!sym(lq), ymax = !!sym(uq)),
      fill = "yellow", colour = "black"
    ) +
    geom_segment(
      aes(x = xmin, xend = xmax, y = Median, yend = Median), colour = "red"
    ) +
    geom_linerange(aes(x = Base, ymin = `10th_Percentile`, ymax = !!sym(lq))) +
    geom_linerange(aes(x = Base, ymin = !!sym(uq), ymax = `90th_Percentile`)) +
    geom_line(aes(Base, Mean, group = Filename), colour = "blue") +
    scale_fill_manual(values = getColours(pwfCols)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    xlab("Position in read (bp)") +
    ylab(ylab) +
    facet_wrap(~Filename, ncol = 1) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

  if (usePlotly) {
    hv <-  c(
      "Base", "Mean", "Median", "Upper_Quartile", "Lower_Quartile",
      "`10th_Percentile`", "`90th_Percentile`"
    )
    qualPlot <- qualPlot +
      xlab("") +
      theme(legend.position = "none")
    qualPlot <- suppressMessages(
      suppressWarnings(
        plotly::ggplotly(qualPlot, hoverinfo = hv)
      )
    )
    qualPlot <- suppressMessages(
      suppressWarnings(
        plotly::subplot(
          plotly::plotly_empty(),
          qualPlot,
          widths = c(0.15,0.85)
        )
      )
    )
    qualPlot <- plotly::layout(qualPlot, yaxis2 = list(title = ylab))
    ## Set the hoverinfo for bg rectangles to none,
    ## This will effectively hide them
    qualPlot$x$data <- lapply(qualPlot$x$data, .hidePWFRects)
    ## Turn off the boxplot fill hover
    qualPlot$x$data[[5]]$hoverinfo <- "none"
    ## Remove xmax & xmin from the hover info
    qualPlot$x$data[[6]]$text <- gsub(
      "xmax:.+(Median.+)xmin.+", "\\1",qualPlot$x$data[[6]]$text
    )

  }

  qualPlot
}
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod("plotBaseQuals", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 25, fail = 20,
    boxWidth = 0.8, plotType = c("heatmap", "boxplot"),
    plotValue = "Mean", cluster = FALSE, dendrogram = FALSE,
    nc = 2, heat_w = 8L, ...){

  ## Get the data
  mod <- "Per_base_sequence_quality"
  df <- getModule(x, mod)
  maxQ <- max(df[["90th_Percentile"]])

  if (!length(df)) {
    msg <- "No Per Base Sequence Quality Module Detected"
    qualPlot <- .emptyPlot(msg)
    if (usePlotly) qualPlot <- ggplotly(qualPlot, tooltip = "")
    return(qualPlot)
  }

  ## Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(.isValidPwf(pwfCols))

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]

  ## Get the Illumina encoding for the axis label
  enc <- getModule(x, "Basic_Statistics")$Encoding[1]
  enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)

  ## Set the plot type & value
  plotType <- match.arg(plotType)
  possVals <- setdiff(colnames(df), c("Filename", "Base"))
  plotValue <- match.arg(plotValue, possVals)
  xlab <- "Position in read (bp)"
  lq <- "Lower_Quartile"
  uq <- "Upper_Quartile"

  if (plotType == "boxplot") {

    ylab <- paste0("Quality Scores (", enc, " encoding)")

    ## Sort out the colours
    pwfCols <- setAlpha(pwfCols, 0.2)

    ## Setup the boxes & the x-axis
    stopifnot(is.numeric(boxWidth))
    df$Base <- factor(df$Base, levels = unique(df$Base))
    df$Position <- as.integer(df$Base)
    df$xmin <- df$Position - boxWidth/2
    df$xmax <- df$Position + boxWidth/2

    ## Set the limits & rectangles
    ylim <- c(0, max(df$`90th_Percentile`) + 1)
    expand_x <- round(0.015*(max(df$Position) - min(df$Position)), 1)
    rects <- tibble(
      xmin = min(df$Position) - expand_x,
      xmax = max(df$Position) + expand_x,
      ymin = c(0, fail, warn),
      ymax = c(fail, warn, max(ylim)),
      Status = c("FAIL", "WARN", "PASS")
    )

    ## Get any theme arguments for dotArgs that have been set
    ## manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    ## Generate the basic plot
    df$Filename <- labels[df$Filename]
    qualPlot <- ggplot(df) +
      geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
      ) +
      geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = !!sym(lq), ymax = !!sym(uq)),
        fill = "yellow", colour = "black"
      ) +
      geom_segment(
        aes(x = xmin, xend = xmax, y = Median, yend = Median), colour = "red"
      ) +
      geom_linerange(
        aes(x = Base, ymin = `10th_Percentile`, ymax = !!sym(lq))
      ) +
      geom_linerange(
        aes(x = Base, ymin = !!sym(uq), ymax = `90th_Percentile`)
      ) +
      geom_line(aes(Base, Mean, group = Filename), colour = "blue") +
      scale_fill_manual(values = getColours(pwfCols)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      xlab(xlab) +
      ylab(ylab) +
      facet_wrap(~Filename, ncol = nc) +
      guides(fill = "none") +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
    if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

    ## Make interactive if required
    if (usePlotly) {
      hv <- c(
        "Base", "Mean", "Median", "Upper_Quartile", "Lower_Quartile",
        "`10th_Percentile`", "`90th_Percentile`"
      )
      qualPlot <- qualPlot + theme(legend.position = "none")
      qualPlot <- suppressMessages(
        plotly::ggplotly(qualPlot, hoverinfo = hv)
      )
      ## Set the hoverinfo for bg rectangles to the vertices
      ## only. This will effectively hide them
      qualPlot$x$data <- lapply(qualPlot$x$data, .hidePWFRects)
      qualPlot$x$data <- lapply(qualPlot$x$data, function(x){
        x$text <- gsub("xmax:.+(Median.+)xmin.+", "\\1", x$text)
        x
      })
    }
    out <- qualPlot
  }

  if (plotType == "heatmap") {

    ## Sort out the data for the Base, allowing for any size bins
    df$Base <- lapply(
      df$Base,
      function(x){
        rng <- as.integer(str_split(x, pattern = "-")[[1]])
        seq(min(rng), max(rng), by = 1L)
      }
    )
    df <- unnest(df, Base)
    df <- df[c("Filename", "Base", plotValue)]
    maxVal <- max(df[[plotValue]], na.rm = TRUE)
    phredMax <- ifelse(maxVal <= warn, max(maxQ, 41), ceiling(maxVal + 1))

    ## Set up the dendrogram & labels
    key <- names(labels)
    clusterDend <- .makeDendro(df, "Filename", "Base", plotValue)
    dx <- ggdendro::dendro_data(clusterDend)
    if (dendrogram | cluster) {
      key <- labels(clusterDend)
    }
    if (!dendrogram) dx$segments <- dx$segments[0,]
    lv <- labels[key]
    df$Filename <- factor(labels[df$Filename], levels = lv)

    cols <- .makePwfGradient(
      vals = na.omit(df[[plotValue]]), pwfCols = pwfCols,
      breaks = c(0, fail, warn, phredMax), passLow = FALSE, na.value = "white"
    )

    ## Start the heatmap
    hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
    qualPlot <- ggplot(
      df, aes(Base, Filename, fill = !!sym(plotValue))
    ) +
      geom_tile() +
      labs(x = xlab, y = c()) +
      ggtitle(stringr::str_to_title(gsub("_", " ", mod))) +
      do.call("scale_fill_gradientn", cols) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expansion(0), position = "right") +
      theme(
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = hj),
        axis.title.x = element_text(hjust = hj),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points")
      )

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
    if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

    ## Add the status sideBar
    status <- subset(getSummary(x), Category == gsub("_", " ", mod))
    status$Filename <- factor(labels[status$Filename], levels = lv)

    out <-
      .prepHeatmap(qualPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

  }

  out

}
)
