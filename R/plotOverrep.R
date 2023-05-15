#' @title Plot a summary of Over-represented Sequences
#'
#' @description Plot a summary of Over-represented Sequences for a set of
#' FASTQC reports
#'
#' @details Percentages are obtained by simply summing those within a report.
#' Any possible double counting by FastQC is ignored for the purposes of a
#' simple approximation.
#'
#' Plots generated from a `FastqcData` object will show the top `n`
#' sequences grouped by their predicted source & coloured by whether the
#' individual sequence would cause a WARN/FAIL.
#'
#' Plots generated from a `FastqcDataList` group sequences by predicted
#' source and summarise as a percentage of the total reads.
#'
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' @param pattern Regex to remove from the end of any filenames
#' @param n The number of sequences to plot from an individual file
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param showPwf Show PASS/WARN/FAIL status on the plot
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param ... Used to pass additional attributes to theme() and between methods
#' @param panel_w Width of main panel on output
#' @param expand.x,expand.y Output from `expansion()` or numeric
#' vectors of length 4. Passed to `scale_*_continuous()`
#' @param paletteName Name of the palette for colouring the possible sources
#' of the overrepresented sequences. Must be a palette name from
#' `RColorBrewer`. Ignored if specifying the scaleFill separately
#' @param scaleFill ggplot scale object
#' @param plotlyLegend Show legend on interactive plots
#'
#'
#' @return A standard ggplot2 object
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
#' # A brief summary across all FastQC reports
#' plotOverrep(fdl)
#'
#' @docType methods
#'
#' @importFrom plotly layout ggplotly
#' @importFrom grDevices rgb
#' @importFrom rlang "!!" sym
#' @import ggplot2
#'
#' @name plotOverrep
#' @rdname plotOverrep-methods
#' @export
setGeneric(
  "plotOverrep",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols, ...
  ){
    standardGeneric("plotOverrep")
  }
)
#' @rdname plotOverrep-methods
#' @export
setMethod(
  "plotOverrep", signature = "ANY",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols, ...
  ){
    .errNotImp(x)
  }
)
#' @rdname plotOverrep-methods
#' @export
setMethod(
  "plotOverrep", signature = "character",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols, ...
  ){
    x <- FastqcDataList(x)
    if (length(x) == 1) x <- x[[1]]
    plotOverrep(x, usePlotly, labels, pattern, pwfCols, ...)
  }
)
#' @rdname plotOverrep-methods
#' @export
setMethod(
  "plotOverrep", signature = "FastqcData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
    n = 10, expand.x = c(0, 0, 0.05, 0), expand.y = c(0, 0.6, 0, 0.6),
    plotlyLegend = FALSE, ...
  ){

    mod <- "Overrepresented_sequences"
    df <- getModule(x, mod)

    if (!length(df)) {
      p <- .emptyPlot("No Overrepresented Sequences")
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(x, labels, pattern)
    labels <- labels[names(labels) %in% df$Filename]
    df$Filename <- labels[df$Filename]

    src <- "Possible_Source"
    df <- dplyr::top_n(df, n, Percentage)
    df$Status <- cut(
      df$Percentage,
      breaks = c(0, 0.1, 1, 100),
      labels = c("PASS", "WARN", "FAIL")
    )
    df[[src]] <- gsub(" \\([0-9]*\\% over [0-9]*bp\\)",  "", df[[src]])
    df$Sequence <- factor(df$Sequence, levels = rev(df$Sequence))
    df$Percentage <- round(df$Percentage, 2) / 100
    df <- droplevels(df)

    ## Check the axis expansion
    stopifnot(is.numeric(expand.x), length(expand.x) == 4)
    stopifnot(is.numeric(expand.y), length(expand.y) == 4)

    ## Set plotting parameters
    xLab <- "Percent of Total Reads (%)"
    yLab <- "Overrepresented Sequence"

    ## Sort out the colours & pass/warn/fail breaks
    if (missing(pwfCols)) pwfCols <- getColours(ngsReports::pwf)
    pwfCols <- pwfCols[names(pwfCols) %in% levels(df$Status)]

    fm <- ". ~ ."
    if (!usePlotly) fm <- "Possible_Source ~ ."

    p <- ggplot(
      df, aes(Sequence, Percentage, fill = Status, label = !!sym(src))
    ) +
      geom_bar(stat = "identity") +
      labs(y = xLab, x = yLab) +
      scale_y_continuous(expand = expand.x, labels = scales::percent) +
      scale_x_discrete(expand = expand.y) +
      facet_grid(as.formula(fm), scales = "free_y", space = "free") +
      theme_bw() +
      coord_flip() +
      scale_fill_manual(values = pwfCols)
    p <- .updateThemeFromDots(p, ...)

    if (usePlotly) {

      ## Add the customisations for plotly
      p <- p + ggtitle(df$Filename[1]) +
        theme(
          axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )
      if (!plotlyLegend) p <- p + theme(legend.position = "none")

      ## Add the empty plot to align in the shiny app
      p <- suppressMessages(suppressWarnings(plotly::ggplotly(p)))
    }

    p

  }
)
#' @rdname plotOverrep-methods
#' @export
setMethod(
  "plotOverrep", signature = "FastqcDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
    showPwf = TRUE, cluster = FALSE, dendrogram = FALSE, scaleFill = NULL,
    paletteName = "Set1", panel_w = 8, expand.x = c(0, 0, 0.05, 0),
    expand.y = rep(0, 4), ...
  ){

    mod <- "Overrepresented_sequences"
    df <- getModule(x, mod)

    if (!length(df)) {
      p <- .emptyPlot("No Overrepresented Sequences")
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(x, labels, pattern)
    labels <- labels[names(labels) %in% df$Filename]

    src <- "Possible_Source"
    df[[src]] <- gsub(" \\([0-9]*\\% over [0-9]*bp\\)", "", df[[src]])
    df <- dplyr::group_by(df, Filename, !!sym(src))
    df <- dplyr::summarise(df, Percentage = sum(Percentage), .groups = "keep")
    df <- dplyr::ungroup(df)
    df$Percentage <- round(df$Percentage, 2) / 100
    lev <- unique(dplyr::arrange(df, Percentage)[[src]])
    df[[src]] <- factor(df[[src]], levels = lev)

    ## Now define the order for a dendrogram if required
    key <- names(labels)
    cols <- c("Filename", src, "Percentage")
    clusterDend <- .makeDendro(df[cols], "Filename", src, "Percentage")
    dx <- ggdendro::dendro_data(clusterDend)
    if (dendrogram | cluster) key <- labels(clusterDend)
    if (!dendrogram) dx$segments <- dx$segments[0,]
    ## Now set everything as factors
    df$Filename <- factor(labels[df$Filename], levels = labels[key])
    maxChar <- max(nchar(levels(df$Filename)))

    ## Check the axis expansion
    stopifnot(is.numeric(expand.x), length(expand.x) == 4)
    stopifnot(is.numeric(expand.y), length(expand.y) == 4)

    ## Define the palette
    if (is.null(scaleFill)) {
      scaleFill <- scale_fill_brewer(palette = "Set1")
      if (paletteName %in% rownames(RColorBrewer::brewer.pal.info)) {
        scaleFill <- scale_fill_brewer(palette = paletteName)
      }
    }
    stopifnot(is(scaleFill, "ScaleDiscrete"))
    stopifnot(scaleFill$aesthetics == "fill")

    xLab <- "Overrepresented Sequences (% of Total)"
    p <- ggplot(
      df, aes(Percentage, Filename, fill = !!sym(src))
    ) +
      geom_bar(stat = "identity") +
      labs(x = xLab,y = c(), fill = src) +
      scale_y_discrete(position = "right", expand = expand.y) +
      scale_x_continuous(expand = expand.x, labels = scales::percent) +
      scaleFill +
      theme_bw()
    if (showPwf)
      p <- p + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))
    p <- .updateThemeFromDots(p, ...)

    ## Prepare the status
    status <- getSummary(x)
    status <- subset(status, Category == "Overrepresented sequences")
    status <- subset(status, Filename %in% key)
    status$Filename <- factor(labels[status$Filename], levels = labels[key])
    if (!showPwf) status <- status[0,]

    if (missing(pwfCols)) pwfCols <- pwf
    .prepHeatmap(p, status, dx$segments, usePlotly, panel_w, pwfCols)

  }
)
