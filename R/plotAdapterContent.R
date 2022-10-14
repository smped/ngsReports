#' @title Draw an Adapter Content Plot
#'
#' @description Draw an Adapter Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the Adapter_Content module from the supplied object and
#' generates a ggplot2 object, with a set of minimal defaults.
#' The output of this function can be further modified using the standard
#' ggplot2 methods.
#'
#' When `x` is a single or FastqcData object line plots will always be
#' drawn for all adapters.
#' Otherwise, users can select line plots or heatmaps.
#' When plotting more than one fastqc file, any undetected adapters will not be
#' shown.
#'
#' An interactive version of the plot can be made by setting `usePlotly`
#' as `TRUE`
#'
#' @param x Can be a `FastqcData`, a `FastqcDataList` or character
#' vector of file paths
#' @param usePlotly `logical`. Output as ggplot2 (default) or plotly
#' object.
#' @param adapterType A regular expression matching the adapter(s) to be
#' plotted. To plot all adapters summed, specify `adapterType = "Total"`.
#' This is the default behaviour.
#' @param plotType `character`. Can only take the values
#' `plotType = "heatmap"` or `plotType = "line"`
#' @param warn,fail The default values for warn and fail are 5 and 10
#' respectively (i.e. percentages)
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param heat_w Width of the heatmap relative to other plot components
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
#' plotAdapterContent(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fqName(fdl))
#' plotAdapterContent(fdl[r1])
#'
#' # Plot just the Universal Adapter
#' # and change the y-axis using ggplot2::scale_y_continuous
#' plotAdapterContent(fdl, adapterType ="Illumina_Universal", plotType = "line") +
#' facet_wrap(~Filename) +
#' guides(colour = "none")
#'
#' @docType methods
#'
#' @import ggplot2
#' @import tibble
#' @import patchwork
#' @importFrom plotly plotly_empty ggplotly
#' @importFrom stats hclust dist
#' @importFrom tidyselect one_of all_of
#' @importFrom dplyr across summarise group_by ungroup
#' @importFrom ggdendro segment
#'
#' @name plotAdapterContent
#' @rdname plotAdapterContent-methods
#' @export
setGeneric("plotAdapterContent", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 10, ...){
  standardGeneric("plotAdapterContent")
}
)
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "ANY", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 10, ...){
  .errNotImp(x)
}
)
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "character", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 10, ...){
  x <- FastqcDataList(x)
  if (length(x) == 1) x <- x[[1]]
  plotAdapterContent(x, usePlotly, labels, pwfCols, warn, fail, ...)
}
)
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 10,  ...){

  df <- getModule(x, "Adapter_Content")

  ## Check for values to plot & return an empty plot if none are found
  valueCols <- setdiff(colnames(df), c("Filename", "Position"))
  msg <- c()
  if (!length(df)) msg <- "No Adapter Content Module Detected"
  if (sum(df[valueCols]) == 0)
    msg <- "No Adapter Content Found in Sequences"
  if (!is.null(msg)) {
    acPlot <- .emptyPlot(msg)
    if (usePlotly) acPlot <- ggplotly(acPlot, tooltip = "")
    return(acPlot)
  }

  ## Set any labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]

  ## Sort out the colours & pass/warn/fail breaks
  if (missing(pwfCols)) pwfCols <- pwf
  stopifnot(.isValidPwf(pwfCols))
  stopifnot(is.numeric(c(warn, fail)))
  stopifnot(all(fail < 100, warn < fail,  warn > 0))

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  ## Change to long form and remove the _ symbols between words
  df <- tidyr::gather(df, key = "Type", value = "Percent", one_of(valueCols))
  df$Type <- gsub("_", " ", df$Type)

  ## Set the positions as a factor
  df$Position <- gsub("([0-9]*)-.+", "\\1", df$Position)
  df$Position <- as.numeric(df$Position)
  ## Round the percent for nicer plotting with plotly
  df$Percent <- round(df$Percent, 2)

  ## Add transparency to background colours & define the rectangles
  pwfCols <- setAlpha(pwfCols, 0.2)
  rng <- structure(range(df$Position), names = c("min", "max"))
  rects <- tibble(
    xmin = 0,
    xmax = rng["max"],
    ymin = c(0, warn, fail),
    ymax = c(warn, fail, 100),
    Status = c("PASS", "WARN", "FAIL")
  )

  ## Create the basic plot
  xLab <- "Position in read (bp)"
  yLab <- "Percent (%)"
  acPlot <- ggplot(df) +
    geom_rect(
      data = rects,
      aes_string(
        xmin = "xmin", xmax = "xmax",
        ymin = "ymin", ymax = "ymax",
        fill = "Status")
    ) +
    geom_line(aes_string(x = "Position", y = "Percent", colour = "Type")) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = getColours(pwfCols)) +
    scale_colour_discrete() +
    facet_wrap(~Filename, ncol = 1) +
    labs(x = xLab, y = yLab) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.title = element_blank()
    )

  ## Add the basic customisations
  if (!is.null(userTheme)) acPlot <- acPlot + userTheme
  ## Make interacive if required
  if (usePlotly) {
    acPlot <- acPlot + theme(legend.position = "none")
    acPlot <- suppressWarnings(
      suppressMessages(
        plotly::subplot(
          plotly::plotly_empty(),
          acPlot,
          widths = c(0.14,0.86)
        )
      )
    )
    acPlot <- plotly::layout(
      acPlot, xaxis2 = list(title = xLab), yaxis2 = list(title = yLab)
    )
    ## Set the hoverinfo for bg rectangles to the vertices only,
    ## This will effectively hide them
    acPlot$x$data <- lapply(acPlot$x$data, .hidePWFRects)
  }

  acPlot

}
)
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 10,
    plotType = c("heatmap", "line"), adapterType = "Total",
    cluster = FALSE, dendrogram = FALSE, heat_w = 8L, ...){

  df <- .tidyAc(x, adapterType)

  ## Check for values to plot & return an empty plot if none are found
  msg <- c()
  if (!length(df))
    msg <- "No Requested Adapter Content Found in Sequences"
  if (!is.null(msg)) {
    acPlot <- .emptyPlot(msg)
    return(acPlot)
  }

  ## Sort out the colours & pass/warn/fail breaks
  if (missing(pwfCols)) pwfCols <- pwf
  stopifnot(.isValidPwf(pwfCols))
  stopifnot(is.numeric(c(warn, fail)))
  stopifnot(all(fail < 100, warn < fail,  warn > 0))
  breaks <- c(0, warn, fail, 100)

  ## Check for valid plotType & labels
  plotType <- match.arg(plotType)
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]

  ## If no adapter content is found fo the selected adapter,
  ## output a simple message. Placing this here handles when combined into
  ## Total_Adapter_Content
  if (max(df$Percent) == 0) {
    msg <- paste("No", unique(df$Type), "found")
    return(.emptyPlot(msg))
  }

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  ## Set the axis label for either plotType
  xLab <- "Position in read (bp)"
  key <- names(labels)

  if (plotType == "heatmap") {

    ## Just use the default order as the key if not clustering
    ## Always turn clustering on if dendrogram = TRUE
    yLab <- c()

    ## Set up the dendrogram
    clusterDend <- .makeDendro(df, "Filename", "Position", "Percent")
    dx <- ggdendro::dendro_data(clusterDend)
    if (dendrogram | cluster) {
      key <- labels(clusterDend)
    }
    if (!dendrogram) dx$segments <- dx$segments[0,]
    df$Filename <- factor(labels[df$Filename], levels = labels[key])

    ## Reset the PWF status using current values & form this plot panel
    ## This needs to be recalculated when using Total AC.
    ## Also make the sideBar
    status <- dplyr::summarise(
      dplyr::group_by(df, Filename), Percent = max(Percent, na.rm = TRUE)
    )
    status$Status <- cut(
      status$Percent, breaks = breaks, include.lowest = TRUE,
      labels = c("PASS", "WARN", "FAIL")
    )

    ## Make the heatmap
    cols <-
      .makePwfGradient(df$Percent, pwf, breaks = breaks, na.value = "white")
    hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
    acPlot <- ggplot(
      df, aes_string("Position", "Filename", fill = "Percent", type = "Type")
    ) +
      geom_tile() +
      ggtitle(unique(df$Type)) +
      labs(x = xLab, y = yLab) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_discrete(expand = c(0, 0), position = "right") +
      do.call("scale_fill_gradientn", cols) +
      theme_bw() +
      theme(
        panel.background = element_blank(),
        plot.title = element_text(hjust = hj),
        axis.title.x = element_text(hjust = hj),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points")
      )

    ## Add custom elements
    if (!is.null(userTheme)) acPlot <- acPlot + userTheme

    out <- .prepHeatmap(acPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

  }

  if (plotType == "line") {

    ## No clustering is required here
    yLab <- "Percent (%)"
    df$Filename <- factor(labels[df$Filename], levels = labels[key])

    ## Set the transparency & position of bg rectangles
    pwfCols <- setAlpha(pwfCols, 0.2)
    rects <- tibble(
      xmin = 0,
      xmax = max(df$Position),
      ymin = c(0, warn, fail),
      ymax = c(warn, fail, 100),
      Status = c("PASS", "WARN", "FAIL")
    )

    ## Create the basic plot
    acPlot <- ggplot(df) +
      geom_rect(
        data = rects,
        aes_string(
          xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax",
          fill = "Status"
        )
      ) +
      geom_line(aes_string("Position", "Percent", colour = "Filename")) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_colour_discrete(labels = labels) +
      scale_fill_manual(values = getColours(pwfCols)) +
      guides(fill = "none") +
      labs(x = xLab, y = yLab) +
      facet_wrap(~Type, ncol = 1) +
      theme_bw()

    if (!is.null(userTheme)) acPlot <- acPlot + userTheme

    ## And draw the plot
    if (usePlotly) {
      acPlot <- acPlot + theme(legend.position = "none")
      acPlot <- suppressMessages(
        plotly::ggplotly(acPlot, hoverinfo = c("x", "y", "colour"))
      )

      ## Set the hoverinfo for bg rectangles to the vertices
      ## only. This will effectively hide them from the mouse
      acPlot$x$data <- lapply(acPlot$x$data, .hidePWFRects)
    }
    out <- acPlot
  }
  out
}
)

#' @title Hide PWF tooltips from line plots
#' @description Hide tooltips from PWF rectangles in line plots
#' @param x plotlyObject$x$data
#' @return plotlyObject$x$data
#' @keywords internal
.hidePWFRects <- function(x){
  ## If there is a name component & it contains
  ## PASS/WARN/FAIL set the hoverinfo to none
  if ("name" %in% names(x)) {
    if (grepl("(PASS|WARN|FAIL)", x$name)) {
      x$hoverinfo <- "none"
    }
  }
  x
}

.tidyAc <- function(x, adapterType = "Total") {

  df <- getModule(x, "Adapter_Content")
  valueCols <- setdiff(colnames(df), c("Filename", "Position"))
  df <- tidyr::pivot_longer(
    df, cols = all_of(valueCols), names_to = "Type", values_to = "Percent"
  )
  adaptOpts <- c("Total_Adapter_Content", valueCols)
  adapterType <- match.arg(adapterType, adaptOpts)

  if (adapterType == "Total_Adapter_Content") {
    ## Sum the adapters by filename& position
    df <- summarise(
      group_by(df, Filename, Position),
      Percent = sum(Percent, na.rm = TRUE), .groups = "drop"
    )
    df <- ungroup(df)
  }
  else{
    df <- dplyr::filter(df, Type == adapterType)
  }
  df$Percent <- round(df$Percent, 2)
  ## Set any binned values to be continuous
  df$Position <- lapply(
    df$Position,
    function(x){
      rng <- as.integer(str_split(x, pattern = "-")[[1]])
      seq(min(rng), max(rng), by = 1L)
    }
  )
  df <- unnest(df, Position)
  df$Type <- gsub("_", " ", adapterType)

  ## Now just keep the three required columns
  df[c("Filename", "Position", "Percent", "Type")]

}

