#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' The raw data from [readTotals()] can otherwise be used to manually
#' create a plot.
#'
#' Duplication levels are based on the value shown on FASTQC reports at the
#' top of the DeDuplicatedTotals plot, which is known to be inaccurate.
#' As it still gives a good guide as to sequence diversity it is included as
#' the default. This can be turned off by setting `duplicated = FALSE`.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' @param pattern Regex used to trim the end of filenames
#' @param duplicated `logical`. Include deduplicated read total estimates
#' to plot charts
#' @param bars If `duplicated = TRUE`, show unique and deduplicated reads
#' as "stacked" or "adjacent".
#' @param barCols Colours for duplicated and unique reads.
#' @param expand.x Output from [`ggplot2::expansion`] controlling x-axis
#' expansion, or a numeric vector of length 4
#' @param scalePaired Scale read totals by 0.5 when paired
#' @param scaleY Scale read totals by this value. The default shows the y-axis
#' in millions
#' @param barFill ScaleDiscrete function to be applied to the plot
#' @param labMin Only show labels for filtering categories higher than this
#' values as a proportion of reads
#' @param labVjust,labNudge Used to place labels within each bar. labNudge is
#' passed to nudge_y within `geom_label()`
#' @param labAlpha,labFill,labSize Passed to `geom_label()` in the alpha, fill
#' and size arguments respectively
#' @param ... Used to pass additional attributes to theme()
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
#' # Plot the Read Totals showing estimated duplicates
#' plotReadTotals(fdl)
#'
#' # Plot the Read Totals without estimated duplicates
#' plotReadTotals(fdl, duplicated = FALSE)
#'
#' @return Returns a ggplot or plotly object
#'
#' @docType methods
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom tidyselect one_of
#' @importFrom forcats fct_inorder
#'
#' @name plotReadTotals
#' @rdname plotReadTotals-methods
#' @export
setGeneric("plotReadTotals", function(x, usePlotly = FALSE, labels, ...){
  standardGeneric("plotReadTotals")
}
)
#' @rdname plotReadTotals-methods
#' @export
setMethod(
  "plotReadTotals", signature = "ANY",
  function(x, usePlotly = FALSE, labels, ...){
    .errNotImp(x)
  })
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels,  pattern = ".(fast|fq|bam|sam).*",
    duplicated = TRUE, bars = c("stacked", "adjacent"),
    barCols = c("red","blue"), expand.x = c(0, 0, 0.02, 0), ...){

  stopifnot(is.logical(duplicated))
  df <- readTotals(x)

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, pattern = pattern, ...)
  df$Filename <- labels[df$Filename]
  df$Filename <- forcats::fct_inorder(df$Filename)

  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  ## Get the colours for the barplot
  barCols <- tryCatch(barCols[seq_len(duplicated + 1)])

  ## Check the axis expansion
  stopifnot(length(expand.x) == 4, is.numeric(expand.x))
  xMax <- max(df$Total_Sequences)
  xLab <- "Read Totals"

  if (!duplicated) {

    rtPlot <- ggplot(df, aes(Filename, Total_Sequences)) +
      geom_bar(stat = "identity", fill = barCols) +
      scale_y_continuous(
        labels = scales::comma, limits = c(0, xMax), expand = expand.x
      ) +
      labs(y = xLab) +
      coord_flip() +
      theme_bw()

    if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme
    if (usePlotly) rtPlot <- plotly::ggplotly(rtPlot)

  }

  if (duplicated) {

    bars <- match.arg(bars)

    ## Add the information to a joined data.frame
    deDup <- getModule(x, "Total_Deduplicated_Percentage")
    deDup$Filename <- labels[deDup$Filename]
    deDup$Filename <- factor(deDup$Filename, levels = levels(df$Filename))
    names(deDup) <- gsub("Total_Deduplicated_", "", names(deDup))
    df <- dplyr::left_join(deDup, df, by = "Filename")

    ##Setup the df for plotting
    types <- c("Unique", "Duplicated")
    df$Unique <- df$Percentage*df$Total_Sequences/100
    df$Unique <- round(df$Unique, 0)
    df$Duplicated <- df$Total_Sequences - df$Unique
    df <- df[c("Filename", types)]
    df <- tidyr::gather(df, "Type", "Total", one_of(types))

    barPos <- c(adjacent = "dodge", stacked = "stack")[bars]

    ## The x-axis expansion needs to be reset for this one
    if (bars == "adjacent") xMax <- max(df$Total)*(1 + expand.x[[1]])

    ## Make the plot
    rtPlot <- ggplot(df, aes(Filename, Total, fill = Type)) +
      geom_bar(stat = "identity", position = barPos) +
      scale_y_continuous(
        labels = scales::comma, limits = c(0, xMax), expand = expand.x
      ) +
      scale_fill_manual(values = barCols) +
      labs(y = xLab) +
      coord_flip() +
      theme_bw()

    ## Add common themes & labels
    if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

    if (usePlotly) {

      # Hide the legend
      rtPlot <- rtPlot + theme(legend.position = "none")
      # Render as a plotly object
      rtPlot <- suppressMessages(
        suppressWarnings(plotly::ggplotly(rtPlot))
      )
    }

  }

  ## Draw the plot
  rtPlot
}
)
#' @importFrom dplyr mutate
#' @importFrom scales percent comma
#' @rdname plotReadTotals-methods
#' @export
setMethod(
  "plotReadTotals", signature = "FastpDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam|sam).*",
    scalePaired  = TRUE, scaleY = 1e6,
    barFill = scale_fill_viridis_d(option = "cividis", direction = -1),
    labMin = 0.05, labVjust = 0.5, labAlpha = 1, labFill = "white",
    labSize = 4, labNudge = 0, ...
  ){

    ## Sort out the summary data
    module <- "Summary"
    data <- getModule(x, module)
    df <- data$Filtering_result
    df <- dplyr::filter(df, total > 0)
    if (!nrow(df)) {
      ## If no filtering was performed, just use the Defore_filtering data
      df <- data$Before_filtering
      stopifnot(nrow(df))
      df$total <- df$total_reads
      df$result <- "total_reads"
      df$rate <- 1
    }

    ## Get the correct order for filling columns
    levels <- sort(rank(df$total))
    levels <- unique(names(levels))
    df$result <- factor(df$result, levels = levels)
    df$result <- forcats::fct_relabel(
      df$result,
      function(x) str_remove_all(str_to_title(gsub("_", " ", x)), " Reads")
    )
    paired <- getModule(x, "paired")
    df <- dplyr::left_join(df, paired, by = "Filename")
    y_col <- ifelse(scaleY == 1, "Reads", paste0("Reads (x", comma(scaleY, 1), ")"))
    df[[y_col]] <- df$total / scaleY
    if (scalePaired) df[[y_col]][df$paired] <- 0.5 * df[[y_col]][df$paired]
    df <- mutate(df, cumsum = cumsum(!!sym(y_col)), .by = Filename)
    df$label_y <- df$cumsum - (1 - labVjust) * df[[y_col]]

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, pattern = pattern, ...)
    df$Filename <- factor(labels[df$Filename], levels = labels)

    ## Add additional columns for a nice plotly figure
    df$Total <- df$total
    if (scalePaired) df$Total[df$paired] <- 0.5 * df$Total[df$paired]
    df$Total <- comma(df$Total)
    df[["% Reads"]] <- percent(df$rate, 0.01)
    names(df) <- gsub("result", "Filtering Result", names(df))
    stopifnot(is(barFill, "ScaleDiscrete"))
    p <- ggplot(
      df,
      aes(
        Filename, !!sym(y_col), fill = !!sym("Filtering Result"),
        total = !!sym("Total"), paired = !!sym("paired"),
        percent = !!sym("% Reads")
      )
    ) +
      geom_col() +
      scale_y_continuous(labels = comma, expand = expansion(c(0, 0.05))) +
      barFill +
      theme_bw()
    p <- .updateThemeFromDots(p, ...)

    if (!usePlotly) {
      p <- p +
        geom_label(
          aes(Filename, y = !!sym("label_y"), label = percent(rate, 0.1)),
          data = dplyr::filter(df, rate >= labMin),
          alpha = labAlpha, fill = labFill, size = labSize, nudge_y = labNudge
        )
    } else {
      p <- p + theme(legend.position = "none")
      tt <- c("Filename", "Filtering Result", "Total", "% Reads")
      p <- plotly::ggplotly(p, tooltip = tt)
    }

    p

  }
)
