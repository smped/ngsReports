#' @title Plot Insert Size Distributions
#'
#' @description
#' Plot the insert size distribution from one of Fastp reports
#'
#' @param x A FastpData or FastpDataList object
#' @param usePlotly `logical`. Generate an interactive plot using plotly
#' @param labels An optional named vector of labels for the file names.
#' All file names must be present in the names of the vector.
#' @param pattern Regex to remove from the end of any filenames
#' @param plotType Determine the plot type. Options vary with the input structure
#' @param counts logical(1) Plot read counts, or percentages (default)
#' @param plotTheme a \link[ggplot2]{theme} object
#' @param expand.x,expand.y Axis expansions
#' @param cluster `logical` default `FALSE`. If set to `TRUE`, data will be
#' clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE` the dendrogram
#' will be displayed.
#' @param heat_w Width of the heatmap relative to other plot components
#' @param scaleFill Continuous scale used to fill heatmap cells. Defaults to the
#' "inferno" palette
#' @param scaleColour Discrete scale for adding line colours
#' @param ... Passed to `geom*` functions during plotting
#'
#' @details
#' Takes a Fastp os a set of Fastp reports and plot insert size distributions.
#' Plots can be drawn as cumulative totals or the default histograms for a
#' single report, and as boxplots or heatmaps for a set of reports
#'
#' @return A ggplot or plotly object
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastp.json.gz", full.names = TRUE)
#' fp <- FastpData(fl)
#' plotInsertSize(
#'   fp, counts = TRUE, fill = "steelblue4",
#'   plotTheme = theme(plot.title = element_text(hjust = 0.5))
#' )
#' plotInsertSize(fp, plotType = "cumulative")
#'
#' @name plotInsertSize
#' @rdname plotInsertSize-methods
#' @export
setGeneric(
  "plotInsertSize",
  function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){
    standardGeneric("plotInsertSize")
  }
)
#' @importFrom tidyr unnest
#' @importFrom rlang sym !!
#' @rdname plotInsertSize-methods
#' @export
setMethod(
  "plotInsertSize", signature = "FastpData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
    plotType = c("histogram", "cumulative"), counts = FALSE,
    plotTheme = theme_get(), expand.x = 0.01, expand.y = c(0, 0.05), ...
  ){

    module <- "Insert_size"
    data <- getModule(x, module)
    df <- unnest(data, "histogram")
    df <- df[c("Filename", "insert_size", "histogram", "freq")]
    maxInsert <- max(dplyr::filter(df, !!sym("freq") > 0)$insert_size)
    unknown <- data[c("Filename", "unknown", "unknown_rate")]
    df <- dplyr::filter(df, !!sym("insert_size") <= maxInsert)
    yval <- c("histogram", "freq")[2 - counts]
    lab_fun <- ifelse(counts, comma, percent)

    labels <- .makeLabels(df, labels, pattern, "Filename")
    df$Filename <- factor(labels[df$Filename], levels = labels)

    plotType <- match.arg(plotType)
    stopifnot(is(plotTheme, "theme"))
    df[["%"]] <- scales::percent(df$freq, 0.01)
    df$Total <- scales::comma(df$histogram)
    df[["Cumulative %"]] <- scales::percent(cumsum(df$freq), 0.01)
    df[["Cumulative Total"]] <- scales::comma(cumsum(df$histogram), 1)
    names(df) <- gsub("insert_size", "Insert Size", names(df))

    if (plotType == "histogram") {

      main <- paste0(unique(df$Filename), ": Insert Size Distribution")
      y_lab <- ifelse(counts, "Total Reads", "Reads (%)")
      p <- ggplot(
        df,
        aes(
          !!sym("Insert Size"), !!sym(yval), name = Filename,
          percent = !!sym("%"), total = Total,
          cumperc = !!sym("Cumulative %"), cumtot = !!sym("Cumulative Total")
        )
      ) +
        geom_col(...) +
        scale_x_continuous(expand = expansion(rep_len(expand.x, 2))) +
        scale_y_continuous(
          labels = lab_fun, expand = expansion(rep_len(expand.y, 2))
        ) +
        ggtitle(main) +
        labs(x = "Insert Size", y = y_lab) +
        theme_bw() +
        plotTheme

    }
    if (plotType == "cumulative") {

      df[[yval]] <- cumsum(df[[yval]])
      y_lab <- ifelse(counts, "Cumulative Reads", " Cumulative Reads (%)")
      main <- paste0(unique(df$Filename), ": Cumulative Insert Size Distribution")
      p <- ggplot(
        df,
        aes(
          !!sym("Insert Size"), !!sym(yval), percent = !!sym("%"), total = Total,
          cumperc = !!sym("Cumulative %"), cumtot = !!sym("Cumulative Total")
        )
      ) +
        geom_line(aes(group = Filename), ...) +
        scale_x_continuous(expand = expansion(rep_len(expand.x, 2))) +
        scale_y_continuous(
          labels = lab_fun, expand = expansion(rep_len(expand.y, 2))
        ) +
        ggtitle(main) +
        labs(x = "Insert Size", y = y_lab) +
        theme_bw() +
        plotTheme

    }

    if (usePlotly) {
      tt <- c("Filename", "x", "Total", "Cumulative Total", "%", "Cumulative %")
      p <- plotly::ggplotly(p, tooltip = tt)
    }

    p

  }
)
#' @importFrom tidyr unnest
#' @importFrom rlang sym !! :=
#' @importFrom dplyr mutate across
#' @importFrom scales percent comma
#' @rdname plotInsertSize-methods
#' @export
setMethod(
  "plotInsertSize", signature = "FastpDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
    plotType = c("heatmap", "line", "cumulative"), plotTheme = theme_get(),
    scaleFill = NULL, scaleColour = NULL, cluster = FALSE, dendrogram = FALSE,
    heat_w = 8, ...
  ){

    module <- "Insert_size"
    data <- getModule(x, module)
    df <- unnest(data, "histogram")
    maxInsert <- max(dplyr::filter(df, !!sym("freq") > 0)$insert_size)
    df <- dplyr::filter(df, !!sym("insert_size") <= maxInsert)

    labels <- .makeLabels(df, labels, pattern, "Filename")
    labels <- labels[names(labels) %in% df$Filename]
    key <- names(labels)

    ## Check Theme & rename df columns for nicer looking plots
    stopifnot(is(plotTheme, "theme"))
    plotType <- match.arg(plotType)
    names(df) <- gsub("insert_size", "Insert Size", names(df))
    names(df) <- gsub("^freq$", "Frequency", names(df))

    ## Needed for either cumulative or line
    if (!is.null(scaleColour)) {
      stopifnot(is(scaleColour, "ScaleDiscrete"))
      stopifnot(scaleColour$aesthetics == "colour")
    }

    if (plotType %in% c("line", "cumulative")) {

      vals <- c(line = "Frequency", cumulative = "Cumulative Frequency")
      yval <- vals[plotType]
      df <- mutate(
        df, !!sym("Cumulative Frequency") := cumsum(!!sym("Frequency")),
        .by = Filename
      )
      df <- mutate(df, across(all_of(unname(vals)), \(x) round(100 * x, 3)))
      df$Filename <- factor(labels[df$Filename], levels = labels[key])
      p <- ggplot(df, aes(!!sym("Insert Size"), !!sym(yval))) +
        geom_line(aes(colour = Filename)) +
        plotTheme +
        scaleColour +
        scale_y_continuous(labels = scales::percent_format(scale = 1))
      p <- .updateThemeFromDots(p, ...)
      if (usePlotly) {
        p <- plotly::ggplotly(p)
      }
    }

    if (plotType == "heatmap") {

      if (is.null(scaleFill)) {
        scaleFill <- scale_fill_viridis_c(labels = percent, option = "inferno")
      }
      stopifnot(is(scaleFill, "ScaleContinuous"))
      stopifnot(scaleFill$aesthetics == "fill")

      n <- length(x)
      if (n > 1) {
        clusterDend <- .makeDendro(df, "Filename", "Insert Size", "Frequency")
        dx <- ggdendro::dendro_data(clusterDend)
        if (dendrogram | cluster) key <- labels(clusterDend)
      } else {
        cluster <- dendrogram <- FALSE
        dx <- list()
        dx$segments <- lapply(rep_len(0, 4), numeric)
        names(dx$segments) <- c("x", "y", "xend", "yend")
        dx$segments <- as_tibble(dx$segments)
      }
      if (!dendrogram) dx$segments <- dx$segments[0,]
      ## Now set everything as factors
      df$Filename <- factor(labels[df$Filename], levels = labels[key])
      df[["%"]] <- percent(df$Frequency, 0.01)
      df$Total <- comma(df$histogram)

      p <- ggplot(
        df,
        aes(
          x = !!sym("Insert Size"), y = Filename, fill = !!sym("Frequency"),
          perc = !!sym("%"), total = !!sym("Total")
        )
      ) +
        geom_raster() +
        ggtitle("Insert Size Distribution") +
        scale_x_continuous(expand = rep_len(0, 4)) +
        scale_y_discrete(expand = rep_len(0, 4), position = "right") +
        scaleFill +
        plotTheme
      if (dendrogram)
        p <- p + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))

      tt <- c("Filename", "Insert Size", "%", "Total")
      p <- .prepHeatmap(
        p, tibble(), usePlotly, heat_w, segments = dx$segments, hv = tt
      )

    }

    p

  }
)
