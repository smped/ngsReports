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
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' @param pattern Regex to remove from the end of the Fastp report and Fastq
#' file names
#' @param nc `numeric`. The number of columns to create in the plot layout.
#' Only used if drawing boxplots for multiple files in a FastqcDataList
#' @param warn,fail The default values for warn and fail are 30 and 20
#' respectively (i.e. percentages)
#' @param plotType `character` Can be either `"boxplot"` or
#' `"heatmap"`
#' @param plotValue `character` Type of data to be presented. Can be
#' any of the columns returned by the appropriate call to `getModule()`
#' @param reads Create plots for read1, read2 or all when using a FastpDataList
#' @param readsBy If paired reads are present, separate using either linetype
#' or by facet
#' @param module Select Before and After filtering when using a FastpDataList
#' @param bases Which bases to include on the plot
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail values in plot
#' @param showPwf Include the Pwf status colours
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param boxWidth set the width of boxes when using a boxplot
#' @param heat_w Relative width of any heatmap plot components
#' @param scaleColour ggplot discrete colour scale, passed to lines
#' @param scaleFill ggplot2 continuous scale. Passed to heatmap cells
#' @param plotTheme \link[ggplot2]{theme} object
#' @param plotlyLegend logical(1) Show legend for interactive plots. Only called
#' when drawing line plots
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
#' # FastpData objects have qyalities by base
#' fp <-  FastpData(system.file("extdata/fastp.json", package = "ngsReports"))
#' plotBaseQuals(
#'   fp, plotTheme = theme(plot.title = element_text(hjust = 0.5))
#' )
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
setGeneric(
  "plotBaseQuals",
  function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){
    standardGeneric("plotBaseQuals")
  }
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod(
  "plotBaseQuals", signature = "ANY",
  function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){
    .errNotImp(x)
  }
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod(
  "plotBaseQuals", signature = "FastqcData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
    warn = 25, fail = 20, boxWidth = 0.8, showPwf = TRUE, plotlyLegend = FALSE,
    ...
  ){

    ## Get the data
    df <- getModule(x, "Per_base_sequence_quality")

    ## Make a blank plot if no data is found
    if (!length(df)) {
      msg <- "No Per Base Sequence Quality Module Detected"
      p <- .emptyPlot(msg)
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    # Get the plot labels organised
    labels <- .makeLabels(x, labels, pattern)
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
    if (!showPwf) rects <- rects[NULL,]

    ## Get the Illumina encoding
    enc <- getModule(x, "Basic_Statistics")$Encoding[1]
    enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)
    ylab <- paste0("Quality Scores (", enc, " encoding)")

    ## Generate the basic plot
    lq <- "Lower_Quartile"
    uq <- "Upper_Quartile"
    p <- ggplot(df)
    if (showPwf) p <- p + geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
    )
    p <- p + geom_rect(
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
    p <- .updateThemeFromDots(p, ...)

    if (usePlotly) {
      hv <-  c(
        "Base", "Mean", "Median", "Upper_Quartile", "Lower_Quartile",
        "`10th_Percentile`", "`90th_Percentile`"
      )
      if (!plotlyLegend) p <- p + theme(legend.position = "none")
      p <- suppressMessages(
        suppressWarnings(plotly::ggplotly(p, tooltip = hv))
      )

      hideRects <- function(x) {
        x_matches <- all(x$x %in% unlist(rects[c("xmax", "xmin")]))
        y_matches <- all(x$y %in% unlist(rects[c("ymin", "ymax")]))
        if (x_matches & y_matches) x$hoverinfo <- "none"
        if (all(x$text == "")) x$hoverinfo <- "none" #Box fill
        x
      }

      ## Set the hoverinfo for bg rectangles to none,
      ## This will effectively hide them
      p$x$data <- lapply(p$x$data, hideRects)

    }

    p
  }
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod(
  "plotBaseQuals", signature = "FastqcDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
    pwfCols, warn = 25, fail = 20, showPwf = TRUE, boxWidth = 0.8,
    plotType = c("heatmap", "boxplot"), plotValue = "Mean",
    cluster = FALSE, dendrogram = FALSE, nc = 2, heat_w = 8L, ...
  ){

    ## Get the data
    mod <- "Per_base_sequence_quality"
    df <- getModule(x, mod)
    maxQ <- max(df[["90th_Percentile"]])

    if (!length(df)) {
      msg <- "No Per Base Sequence Quality Module Detected"
      p <- .emptyPlot(msg)
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    ## Sort out the colours
    if (missing(pwfCols)) pwfCols <- ngsReports::pwf
    stopifnot(.isValidPwf(pwfCols))

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(x, labels, pattern)
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
        xmin = min(df$Position) - expand_x, xmax = max(df$Position) + expand_x,
        ymin = c(0, fail, warn), ymax = c(fail, warn, max(ylim)),
        Status = c("FAIL", "WARN", "PASS")
      )

      ## Generate the basic plot
      df$Filename <- labels[df$Filename]
      p <- ggplot(df)
      if (showPwf) p <- p + geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status)
      )
      p <- p + geom_rect(
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
      p <- .updateThemeFromDots(p, ...)

      ## Make interactive if required
      if (usePlotly) {
        hv <- c(
          "Base", "Mean", "Median", "Upper_Quartile", "Lower_Quartile",
          "`10th_Percentile`", "`90th_Percentile`"
        )
        p <- p + theme(legend.position = "none")
        p <- suppressMessages(plotly::ggplotly(p, tooltip = hv))
        hideRects <- function(x) {
          x_matches <- all(x$x %in% unlist(rects[c("xmax", "xmin")]))
          y_matches <- all(x$y %in% unlist(rects[c("ymin", "ymax")]))
          if (x_matches & y_matches) x$hoverinfo <- "none"
          if (all(x$text == "")) x$hoverinfo <- "none" #Box fill
          x
        }
        p$x$data <- lapply(p$x$data, hideRects)

      }
      out <- p
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
      p <- ggplot(df, aes(Base, Filename, fill = !!sym(plotValue))) +
        geom_tile() +
        labs(x = xlab, y = c()) +
        ggtitle(stringr::str_to_title(gsub("_", " ", mod))) +
        do.call("scale_fill_gradientn", cols) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(expand = expansion(0), position = "right") +
        theme(
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = hj),
          axis.title.x = element_text(hjust = hj)
        )
      if (dendrogram | showPwf)
        p <- p + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))

      ## Get any arguments for dotArgs that have been set manually
      p <- .updateThemeFromDots(p, ...)

      ## Add the status sideBar
      status <- subset(getSummary(x), Category == gsub("_", " ", mod))
      status$Filename <- factor(labels[status$Filename], levels = lv)
      if (!showPwf) status <- status[0, ]

      out <- .prepHeatmap(p, status, dx$segments, usePlotly, heat_w, pwfCols)

    }

    out

  }
)
#' @importFrom tidyselect everything all_of
#' @importFrom tidyr pivot_longer
#' @importFrom rlang "!!" sym
#' @rdname plotBaseQuals-methods
#' @export
setMethod(
  "plotBaseQuals", signature = "FastpData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
    pwfCols, warn = 25, fail = 20, showPwf = FALSE,
    module = c("Before_filtering", "After_filtering"),
    reads = c("read1", "read2"), readsBy = c("facet", "linetype"),
    bases = c("A", "T", "C", "G", "mean"),
    scaleColour = NULL, plotTheme = theme(), plotlyLegend = FALSE, ...
  ){

    ## Get the data
    module <- match.arg(module)
    data <- getModule(x, module)
    avail_reads <- names(data)
    reads <- match.arg(reads, avail_reads, several.ok = TRUE)
    df <- dplyr::bind_rows(data[reads], .id = "reads")
    fpName <- .makeLabels(df, labels, pattern)
    df[["Filename"]] <- df[["fqName"]]
    df <- df[c("Filename", "reads", "quality_curves")]
    df <- tidyr::unnest(df, everything())
    bases <- match.arg(bases, colnames(df), several.ok = TRUE)
    readsBy <- match.arg(readsBy)
    linetype <- NULL
    if (readsBy == "linetype") linetype <- sym("Filename")

    ## Make a blank plot if no data is found
    if (!length(df)) {
      msg <- "No Quality Curve Data Detected"
      p <- .emptyPlot(msg)
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    # Get the plot labels organised
    labels <- .makeLabels(df, labels, pattern = pattern, ...)
    labels <- labels[names(labels) %in% df$Filename]
    df$Filename <- labels[df$Filename]
    df <- pivot_longer(
      df, cols = all_of(bases), names_to = "base", values_to = "quality"
    )
    df[["base"]] <- factor(df[["base"]], levels = bases)
    df[["quality"]] <- round(df[["quality"]], 2)

    ## Sort out the Pwf Colours
    if (missing(pwfCols)) pwfCols <- ngsReports::pwf
    stopifnot(.isValidPwf(pwfCols))
    pwfCols <- setAlpha(pwfCols, 0.2)

    ## Set the limits & rectangles
    ylim <- c(0.975 * min(df$quality), 1.025 * max(df$quality))
    expand_x <- round(0.015*(max(df$position) - min(df$position)), 1)
    rects <- tibble(
      xmin = min(df$position) - expand_x,
      xmax = max(df$position) + expand_x,
      ymin = c(0, fail, warn),
      ymax = c(fail, warn, max(ylim)),
      Status = c("FAIL", "WARN", "PASS")
    )
    xlim <- c(rects$xmin[[1]], rects$xmax[[1]])
    rects <- dplyr::filter(rects, !!sym("ymax") > min(ylim))
    rects[["ymin"]][rects[["ymin"]] < min(ylim)] <- min(ylim)

    ## The scaleColour
    if (is.null(scaleColour)) {
      line_cols = c(
        A = "#807C58", `T` = "#601490", C = "green", G = "blue", mean = "grey30"
      )[bases]
      scaleColour <- scale_colour_manual(values = line_cols)
    }
    stopifnot(is(scaleColour, "ScaleDiscrete"))
    stopifnot(scaleColour$aesthetics == "colour")
    stopifnot(is(plotTheme, "theme"))

    names(df) <- stringr::str_to_title(names(df))
    p <- ggplot(df)
    if (showPwf) {
      p <- p + geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Status),
        inherit.aes = FALSE
      )
    }

    p <- p +
      geom_line(
        aes(
          !!sym("Position"), !!sym("Quality"), colour = !!sym("Base"),
          linetype = {{ linetype }}
        ), ...
      ) +
      scaleColour +
      scale_fill_manual(values = getColours(pwfCols)) +
      scale_x_continuous(expand = rep(0, 4), limits = xlim) +
      scale_y_continuous(expand = rep(0, 4), limits = ylim) +
      labs(colour = c()) +
      guides(fill = "none") +
      ggtitle(fpName) +
      theme_bw() +
      plotTheme
    if (readsBy == "facet") p <- p + facet_wrap(~Filename)

    if (usePlotly) {
      tt <- c("Position", "Base", "Quality")
      if (!plotlyLegend) p <- p + theme(legend.position = "none")
      p <- suppressMessages(
        suppressWarnings(plotly::ggplotly(p, tooltip = tt))
      )
      p$x$data <- lapply(
        p$x$data,
        function(x) {
          if (all(x$text == "")) x$hoverinfo <- "none"
          x
        }
      )

    }

    p
  }
)
#' @rdname plotBaseQuals-methods
#' @export
setMethod(
  "plotBaseQuals", signature = "FastpDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
    pwfCols, warn = 25, fail = 20, showPwf = FALSE,
    module = c("Before_filtering", "After_filtering"),
    plotType = "heatmap", plotValue = c("mean", "A", "T", "C", "G"),
    scaleFill = NULL, plotTheme = theme(), cluster = FALSE,
    dendrogram = FALSE, heat_w = 8L, ...
  ){

    ## Get the data
    module <- match.arg(module)
    plotType <- match.arg(plotType)
    plotValue <- match.arg(plotValue)

    df <- bind_rows(getModule(x, module), .id = "reads")
    df$Filename <- df$fqName
    df <- df[c("Filename", "quality_curves")]
    df <- unnest(df, !!sym("quality_curves"))
    df <- df[c("Filename", "position", plotValue)]
    maxQ <- max(df[[plotValue]], na.rm = TRUE)
    phredMax <- ifelse(maxQ <= warn, max(maxQ, 41), ceiling(maxQ + 1))

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, pattern)
    labels <- labels[names(labels) %in% df$Filename]
    ## Set up the dendrogram & labels
    key <- names(labels)
    clusterDend <- .makeDendro(df, "Filename", "position", plotValue)
    dx <- ggdendro::dendro_data(clusterDend)
    if (dendrogram | cluster) key <- labels(clusterDend)
    if (!dendrogram) dx$segments <- dx$segments[0,]
    lv <- labels[key]
    df$Filename <- factor(labels[df$Filename], levels = lv)
    df[[plotValue]] <- round(df[[plotValue]], 2)

    ## Set the plot titles
    xlab <- "Position in read (bp)"
    main <- gsub("_", " ", module)
    main <- stringr::str_to_title(paste(main, collapse = ""))

    ## Sort out the colours
    if (missing(pwfCols)) pwfCols <- ngsReports::pwf
    stopifnot(.isValidPwf(pwfCols))
    if (is.null(scaleFill)) {
      cols <- .makePwfGradient(
        vals = na.omit(df[[plotValue]]), pwfCols = pwfCols,
        breaks = c(0, fail, warn, phredMax), passLow = FALSE, na.value = "white"
      )
      scaleFill <- do.call("scale_fill_gradientn", cols)
    }
    stopifnot(is(scaleFill, "ScaleContinuous"))
    stopifnot(scaleFill$aesthetics == "fill")
    stopifnot(is(plotTheme, "theme"))

    ## Start the heatmap
    hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
    p <- ggplot(
      df, aes(!!sym("position"), Filename, fill = !!sym(plotValue))
    ) +
      geom_tile() +
      labs(x = xlab, y = c()) +
      ggtitle(main) +
      scaleFill +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expansion(0), position = "right") +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(), panel.background = element_blank()
      ) +
      plotTheme
    if (dendrogram) p <- p +
      theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))

    ## Manually replicate the status using the provided values
    status <- lapply(split(df, df$Filename), function(x) min(x[[plotValue]]))
    status <- unlist(status)
    status_df <- tibble(Filename = names(status), value = as.numeric(status))
    status_df$Status <- cut(
      status_df$value, breaks = c(0, fail, warn, Inf),
      labels = c("FAIL", "WARN", "PASS")
    )
    status_df$Filename <- factor(status_df$Filename, levels = lv)
    if (!showPwf) status_df <- status_df[0, ]
    .prepHeatmap(p, status_df, dx$segments, usePlotly, heat_w, pwfCols)

  }
)

