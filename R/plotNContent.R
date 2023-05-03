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
#' @param pattern Regex used to trim filenames
#' @param module Used for Fastp* structures to show results before or after
#' filtering
#' @param reads Show plots for read1, read2 or both.
#' @param lineCol Defaults to red
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param heat_w Relative width of any heatmap plot components
#' @param fill continuous scale for ggplot objects
#' @param lineCol Line colours
#' @param colour,linetype Data arguments to decorate plots. Should be drawn
#' from columns in the formatted data, such as "module", "reads", "fqName"
#' @param facetBy Formula to facet the plot by. Passed to `facet_wrap()`
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
    x, usePlotly = FALSE, labels, ...){
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
    xmin = 0, xmax = max(df$xValue),
    ymin = c(0, warn, fail), ymax = c(warn, fail, 100),
    Status = c("PASS", "WARN", "FAIL")
  )

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
  nPlot <- .updateThemeFromDots(nPlot, ...)

  if (usePlotly) {
    nPlot <- nPlot +
      xlab("") +
      theme(legend.position = "none")
    nPlot <- suppressMessages(plotly::ggplotly(nPlot))
    # nPlot <- suppressMessages(
    #   suppressWarnings(
    #     plotly::subplot(
    #       plotly::plotly_empty(),
    #       nPlot,
    #       widths = c(0.14,0.86)
    #     )
    #   )
    # )
    nPlot <- plotly::layout(nPlot, yaxis1 = list(title = yLab))


    ## Set the hoverinfo for bg rectangles to the vertices only,
    ## This will effectively hide them
    nPlot$x$data <- lapply(nPlot$x$data, .hidePWFRects)
    ## Hide the xValue parameter to make it look nicer
    nPlot$x$data[[6]]$text <- gsub(
      "(.+)(xValue.+)(Percentage.+)", "\\1\\3", nPlot$x$data[[6]]$text
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
    breaks = c(0, warn, fail, 101), passLow = TRUE, na.value = "white"
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
  nPlot <- .updateThemeFromDots(nPlot, ...)

  ## Reset the status using current values
  status <- dplyr::summarise(
    dplyr::group_by(df, Filename), Percentage = max(Percentage, na.rm = TRUE)
  )
  status$Status <- cut(
    status$Percentage, include.lowest = TRUE,
    breaks = c(0, warn, fail, 101), labels = c("PASS", "WARN", "FAIL")
  )
  status <- subset(status, Filename %in% key)
  status$Filename <- factor(labels[status$Filename], levels = labels[key])

  .prepHeatmap(nPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

}
)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr unnest
#' @importFrom tidyselect any_of
#' @importFrom scales label_percent comma
#' @importFrom rlang sym "!!"
#' @rdname plotNContent-methods
#' @export
setMethod(
  "plotNContent", signature = "FastpData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".fastq.*|.fastp.*|.fq*",
    module = c("Before_filtering", "After_filtering"),
    reads = c("read1", "read2"), linetype = NULL, colour = NULL,
    lineCol = c("navyblue", "red3"), facetBy = module ~ fqName, ...
  ){

    ## Check args
    mod <- match.arg(module, several.ok = TRUE) # We can plot B4/After
    reads <- match.arg(reads, several.ok = TRUE)

    data <- lapply(mod, function(m) getModule(x, m))
    names(data) <- mod
    data <- lapply(data, bind_rows, .id = "reads")
    df <- bind_rows(data,.id = "Module")
    cols <- c(
      "Module", "reads", "Filename", "fqName", "total_reads", "content_curves",
      "N", "position"
    )
    df <- dplyr::select(df, any_of(cols))
    df <- unnest(df, !!sym("content_curves"))
    df <- dplyr::select(df, any_of(cols))

    ## Drop the suffix, or check the alternate labels
    lbl_df <- dplyr::distinct(df, Filename, fqName)
    fnames <- .makeLabels(lbl_df, labels, pattern = pattern, ...)
    df$Filename <- factor(fnames[df$Filename], levels = fnames)
    df$Filename <- fct_rev(df$Filename)
    fqNames <- .makeLabels(lbl_df, labels, pattern = pattern, col = "fqName")
    df$fqName <- fqNames[df$fqName]
    df$Module <- factor(df$Module, levels = mod)
    df$module <- df$Module ## Makes it easier for users to make typos

    ## Make a blank plot if no data is found
    if (!length(df)) {
      msg <- "No N Content in Reports"
      p <- .emptyPlot(msg)
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    df[["% N"]] <- round(100 * df$N, 2)
    df[["N Reads"]] <- comma(df$N * df$total_reads, 1)
    names(df) <- gsub("position", "Position", names(df))
    df_cols <- colnames(df)
    if (!is.null(linetype)) linetype <- sym(match.arg(linetype, df_cols))
    if (!is.null(colour)) colour <- sym(match.arg(colour, df_cols))
    p <- ggplot(df) +
      geom_line(
        aes(
          Position, !!sym("% N"), colour = {{ colour }}, linetype = {{ linetype }}
        )
      ) +
      facet_grid(facetBy) +
      scale_colour_manual(values = lineCol) +
      scale_y_continuous(
        labels = label_percent(scale = 1), expand = expansion(c(0.01, 0.05))
      ) +
      labs(y = "% N") +
      theme_bw()
    p <- .updateThemeFromDots(p, ...)

    if (usePlotly) {
      ############################################
      ## This is problematic & needs work still ##
      ############################################
      p <- p + theme(legend.position = "none")
      hv <- c("fqName", "Position", "% N", "N Reads", "colour")
      p <- suppressWarnings(plotly::ggplotly(p, tooltip = hv))
    }
    p
  }
)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr unnest
#' @importFrom tidyselect any_of
#' @importFrom scales percent comma
#' @importFrom rlang sym "!!"
#' @rdname plotNContent-methods
#' @export
setMethod(
  "plotNContent", signature = "FastpDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".fastq.*|.fastp.*|.fq*",
    module = c("Before_filtering", "After_filtering"),
    reads = c("read1", "read2"), fill = scale_fill_viridis_c(), ...
  ){


    ## Check args
    mod <- match.arg(module, several.ok = TRUE) # We can plot B4/After
    reads <- match.arg(reads, several.ok = TRUE)
    stopifnot(is(fill, "ScaleContinuous"))

    ## Setup the data
    data <- lapply(mod, function(m) getModule(x, m)[reads])
    names(data) <- mod
    data <- lapply(data, bind_rows, .id = "reads")
    df <- bind_rows(data,.id = "Module")
    cols <- c(
      "Module", "reads", "Filename", "fqName", "total_reads", "content_curves",
      "N", "position"
    )
    df <- dplyr::select(df, any_of(cols))
    df <- unnest(df, !!sym("content_curves"))
    df <- dplyr::select(df, any_of(cols))

    ## Drop the suffix, or check the alternate labels
    lbl_df <- dplyr::distinct(df, Filename, fqName)
    fnames <- .makeLabels(lbl_df, labels, pattern = pattern, ...)
    df$Filename <- factor(fnames[df$Filename], levels = fnames)
    df$Filename <- forcats::fct_rev(df$Filename)
    fqNames <- .makeLabels(lbl_df, labels, pattern = pattern, col = "fqName")
    df$fqName <- fqNames[df$fqName]

    ## Make a blank plot if no data is found
    if (!length(df)) {
      msg <- "No N Content in Reports"
      p <- .emptyPlot(msg)
      if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
      return(p)
    }

    ## Tidy up for plotting
    df$Module <- factor(df$Module, levels = mod)
    df[["% N"]] <- percent(df$N, 0.01)
    df[["N Reads"]] <- comma(df$N * df$total_reads, 1)
    names(df) <- gsub("position", "Position", names(df))
    fm <- as.formula("Module ~ reads")
    p <- ggplot(
      df,
      aes(
        Position, Filename, fill = !!sym("N"), label = !!sym("fqName"),
        percent = !!sym("% N"), total = !!sym("N Reads")
      )
    ) +
      geom_raster() +
      facet_grid(fm) +
      scale_fill_viridis_c(labels = percent) +
      scale_y_discrete(expand = rep(0, 4)) +
      scale_x_continuous(expand = rep(0, 4))
    p <- .updateThemeFromDots(p, ...)

    if (usePlotly) {
      p <- p + theme(legend.position = "none")
      hv <- c("fqName", "Position", "% N", "N Reads", "Module")
      p <- suppressWarnings(plotly::ggplotly(p, tooltip = hv))
    }
    p
  }
)
