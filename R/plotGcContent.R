#' @title Plot the Per Sequence GC Content
#'
#' @description Plot the Per Sequence GC Content for a set of FASTQC files
#'
#' @details
#' Makes plots for GC_Content.
#' When applied to a single FastqcData object a simple line plot will be drawn,
#' with Theoretical GC content overlaid if desired.
#'
#' When applied to multiple FastQC reports, the density at each GC content bin
#' can be shown as a heatmap by setting `theoreticalGC = FALSE`. By
#' default the difference in observed and expected theoretical GC is shown.
#' Species and genome/transcriptome should also be set if utilising the
#' theoretical GC content.
#'
#' As an alternative to a heatmap, a series of overlaid distributions can be
#' shown by setting `plotType = "line"`.
#'
#' Can produce a static ggplot2 object or an interactive plotly object.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or character
#' vector of file paths
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param counts `logical`. Plot the counts from each file if
#' `counts = TRUE`, otherwise frequencies will be plotted.
#' Ignored if calling the function on a FastqcDataList.
#' @param theoreticalGC `logical` default is `FALSE` to give the true
#' GC content, set to `TRUE` to normalize values of GC_Content by the
#' theoretical values using [gcTheoretical()]. `species` must be
#' specified. For Fastqc* objects, the entire distributions will be used,
#' wheras for the Fastp* objects, only the expected mean value is shown as a
#' horizontal line
#' @param gcType `character` Select type of data to normalize GC
#' content against. Accepts either "Genome" (default) or "Transcriptome".
#' @param GCobject an object of class GCTheoretical.
#'  Defaults to the gcTheoretical object supplied with the package
#' @param Fastafile a fasta file contains DNA sequences to generate theoretical
#' GC content
#' @param n number of simulated reads to generate theoretical GC content from
#' `Fastafile`
#' @param species `character` if `gcTheory` is `TRUE` it must be
#' accompanied by a species. Species currently supported can be obtained using
#' `mData(gcTheoretical)`
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param pattern Pattern to remove from the end of filenames
#' @param plotType Takes values "line", "heatmap" or "cdf"
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail values in plot
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram  will be displayed.
#' @param lineCol Colors for observed and theoretical GC lines in single plots
#' @param fillCol Bar colours when plotting data from a FastpDataList
#' @param heat_w Relative width of any heatmap plot components
#' @param ... Used to pass various potting parameters to theme.
#'
#' @return A ggplot2 or plotly object
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
#' # The default plot for a FastqcDataList
#' plotGcContent(fdl)
#'
#' # Plot a single FastqcData object
#' plotGcContent(fdl[[1]])
#'
#' @docType methods
#'
#' @importFrom grDevices colorRampPalette hcl.colors
#' @importFrom stats hclust dist
#' @import ggplot2
#' @importFrom stringr str_to_title
#' @importFrom rlang "!!" sym
#' @name plotGcContent
#' @rdname plotGcContent-methods
#' @export
setGeneric("plotGcContent", function(x, usePlotly = FALSE, labels, ...){
  standardGeneric("plotGcContent")
}
)
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "ANY", function(
    x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
    gcType = c("Genome", "Transcriptome"), species = "Hsapiens",
    GCobject, Fastafile, n = 1e+6, ...){
  .errNotImp(x)
}
)
#' @rdname plotGcContent-methods
#' @export
setMethod(
  "plotGcContent", signature = "FastqcData",
  function(
    x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
    gcType = c("Genome", "Transcriptome"), species = "Hsapiens", GCobject,
    Fastafile, n = 1e6, counts = FALSE, lineCol = c("red3", "black"),
    ...
  ){

  df <- getModule(x, "Per_sequence_GC_content")

  if (!length(df)) {
    gcPlot <- .emptyPlot("No Duplication Levels Module Detected")
    if (usePlotly) gcPlot <- ggplotly(gcPlot, tooltip = "")
    return(gcPlot)
  }

  df$Type <- "GC count per read"

  ## Get the correct y-axis label
  yLab <- c("Frequency", "Count")[counts + 1]

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]

  ## Tidy up the GC content variables
  if (missing(GCobject)) GCobject <- ngsReports::gcTheoretical

  gcTheoryDF <- c()
  subTitle <- c()
  if (theoreticalGC) {
    if (!missing(Fastafile)) {
      rdLen <- max(getModule(x, "Basic_Statistics")$Longest_sequence)
      gcTheoryDF <- estGcDistn(Fastafile, n, rdLen)
      subTitle <-
        paste("Theoretical Distribution based on file", Fastafile)
    }
    else{
      gcType <- stringr::str_to_title(gcType)
      gcType <- match.arg(gcType)
      avail <- gcAvail(GCobject, gcType)
      species <- match.arg(species, avail$Name)
      gcTheoryDF <- getGC(GCobject, species, gcType)
      names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
      subTitle <- paste(
        "Theoretical Distribution based on the",
        species,
        gcType
      )
    }
    gcTheoryDF$Type <- "Theoretical Distribution"
    gcTheoryDF$Filename <- "Theoretical Distribution"
    gcTheoryDF$Freq <- round(gcTheoryDF$Freq,4)
  }

  xLab <- "GC Content (%)"

  if (!counts) {## If using frequencies (not counts)

    vals <- c("GC_Content", "Freq", "Type")
    ## Summarise to frequencies & initialise the plot
    df$Freq <- df$Count / sum(df$Count)
    df <- df[vals]
    df <- dplyr::bind_rows(df, gcTheoryDF)
    df$Type <- as.factor(df$Type)
    df$Freq <- round(df$Freq, 4)

    gcPlot <- ggplot(
      df, aes(!!sym(vals[[1]]), !!sym(vals[[2]]), colour = !!sym(vals[[3]]))
    ) +
      geom_line()
  }
  else{

    vals <- c("GC_Content", "Count", "Type")
    df <- df[vals]
    if (theoreticalGC) {
      gcTheoryDF$Count <- gcTheoryDF$Freq * sum(df$Count)
      gcTheoryDF <- gcTheoryDF[vals]
      df <- dplyr::bind_rows(df, gcTheoryDF)
    }
    ## Initialise the plot using counts
    gcPlot <- ggplot(
      df, aes(!!sym(vals[[1]]), !!sym(vals[[2]]), colour = !!sym(vals[[3]]))
    ) +
      geom_line()
  }

  gcPlot <- gcPlot +
    scale_colour_manual(values = lineCol) +
    scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0.02, 0)) +
    labs(x = xLab, y = yLab, colour = c()) +
    ggtitle(label = labels, subtitle = subTitle) +
    theme_bw() +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_rect(colour = "grey20", linewidth = 0.2),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  gcPlot <- .updateThemeFromDots(gcPlot, ...)

  if (usePlotly) {
    value <- c("Freq", "Count")[counts + 1]
    ttip <- c("GC_Content", value, "Type")
    gcPlot <- gcPlot +
      ggtitle(label = labels, subtitle = c()) +
      theme(legend.position = "none")
    gcPlot <- suppressWarnings(
      suppressMessages(plotly::ggplotly(gcPlot, tooltip = ttip)
      )
    )
    # gcPlot <- suppressMessages(
    #   plotly::subplot(
    #     plotly::plotly_empty(), gcPlot, widths = c(0.14,0.86)
    #   )
    # )
    gcPlot <- plotly::layout(
      gcPlot,
      xaxis1 = list(title = xLab),
      yaxis1 = list(title = yLab)
    )

  }

  ## Draw the plot
  gcPlot
}
)
#' @rdname plotGcContent-methods
#' @export
setMethod("plotGcContent", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, theoreticalGC = TRUE,
    gcType = c("Genome", "Transcriptome"), species = "Hsapiens",
    GCobject, Fastafile, n=1e+6, plotType = c("heatmap", "line", "cdf"),
    pwfCols, cluster = FALSE, dendrogram = FALSE, heat_w = 8, ...){

  mod <- "Per_sequence_GC_content"
  df <- getModule(x, mod)

  if (!length(df)) {
    gcPlot <- .emptyPlot("No Duplication Levels Module Detected")
    if (usePlotly) gcPlot <- ggplotly(gcPlot, tooltip = "")
    return(gcPlot)
  }

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, ...)
  labels <- labels[names(labels) %in% df$Filename]

  ## Always use frequencies not counts
  df <- lapply(split(df, f = df$Filename), function(x){
    x$Percent <- 100*x$Count / sum(x$Count)
    x
  })
  df <- dplyr::bind_rows(df)
  df <- df[c("Filename", "GC_Content", "Percent")]

  ## Initialise objects then fill if required
  gcTheoryDF <- c()
  subTitle <- c()
  if (theoreticalGC) {
    if (!missing(Fastafile)) {
      rdLen <- max(getModule(x, "Basic_Statistics")$Longest_sequence)
      gcTheoryDF <- estGcDistn(Fastafile, n, rdLen)
      subTitle <-
        paste("Theoretical Distribution based on", basename(Fastafile))
    }
    else {
      ## Tidy up the GC content variables
      if (missing(GCobject)) GCobject <- gcTheoretical
      gcType <- stringr::str_to_title(gcType)
      gcType <- match.arg(gcType)
      avail <- gcAvail(GCobject, gcType)
      species <- match.arg(species, avail$Name)
      gcTheoryDF <- getGC(
        GCobject, name = species, type = gcType
      )
      names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"
      subTitle <- paste(
        "Theoretical Distribution based on the",
        species, gcType
      )
    }
    gcTheoryDF$Filename <- "Theoretical Distribution"
    gcTheoryDF$Percent <- round(100*gcTheoryDF$Freq,4)
  }

  ## Check for valid plotType arguments & set the x-label
  plotType <- match.arg(plotType)
  xLab <- "GC Content (%)"

  if (plotType %in% c("cdf", "line")) {

    ## Setup a palette with black as the first colour.
    ## Use the paired palette for easier visualisation of paired data
    ## Onlt used for plotType = "line" or plotType = "cdf
    n <- length(x)
    lineCol <- RColorBrewer::brewer.pal(min(12, n), "Paired")
    lineCol <- colorRampPalette(lineCol)(n)
    lineCol <- c("#000000", lineCol)

    ## Select axis labels and plotting values
    ylab <- c(cdf = "Cumulative (%)", line = "Reads (%)")[plotType]

    ## Set the Filename labels & add the TheoreticalGC df
    df$Filename <- labels[df$Filename]
    df <- dplyr::bind_rows(gcTheoryDF, df)
    df$Filename <- factor(df$Filename, levels = unique(df$Filename))

    ## Calculate the CDF then rejoin if required
    if (plotType == "cdf") {
      df <- lapply(
        split(df, f = df$Filename),
        function(y){
          y[["Percent"]] <- cumsum(y[["Percent"]])
          y
        }
      )
      df <- dplyr::bind_rows(df)
    }

    ## Round down to 2 decimal places for playing more nicely with plotly
    df[["Percent"]] <- round(df[["Percent"]], 2)

    gcPlot <- ggplot(
      df, aes(!!sym("GC_Content"), Percent, colour = Filename)
    ) +
      geom_line() +
      scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0.02, 0)) +
      scale_y_continuous(labels = .addPercent) +
      scale_colour_manual(values = lineCol) +
      labs(x = xLab, y = ylab, colour = c()) +
      ggtitle(label = c(), subtitle = subTitle) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    gcPlot <- .updateThemeFromDots(gcPlot, ...)

    if (usePlotly) {

      gcPlot <- gcPlot +
        labs(colour = "Filename") +
        theme(legend.position = "none") +
        ggtitle(NULL)
      gcPlot <- suppressWarnings(
        suppressMessages(
          plotly::ggplotly(
            gcPlot,
            tooltip = c("x", "y", "colour")
          )
        )
      )
    }

  }

  if (plotType == "heatmap") {

    fillLab <- "Frequency" # Default label
    if (theoreticalGC) {
      ## If using theoretical GC, just show the difference
      df <- lapply(split(df, df$Filename), function(x){
        x$Percent <- x$Percent - unlist(gcTheoryDF$Percent)
        x$Percent <- round(x$Percent, 2)
        x
      })
      df <- dplyr::bind_rows(df)
      fillLab <- "Difference"
    }

    ttl <- ifelse(
      theoreticalGC, "Difference from Theoretical GC", gsub("_", " ", mod)
    )
    key <- names(labels)
    clusterDend <- .makeDendro(df, "Filename", "GC_Content", "Percent")
    dx <- ggdendro::dendro_data(clusterDend)
    if (dendrogram | cluster) key <- labels(clusterDend)
    ## Now set everything as factors
    df$Filename <- factor(labels[df$Filename], levels = labels[key])
    if (!dendrogram) dx$segments <- dx$segments[0,]

    ## Draw the heatmap
    hj <- 0.5 * heat_w / (heat_w + 1)
    gcPlot <- ggplot(df, aes(!!sym("GC_Content"), Filename, fill = Percent)) +
      geom_tile() +
      ggtitle(ttl) +
      labs(x = xLab, fill = fillLab, y = c()) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0), position = "right") +
      scale_fill_gradient2(
        low = "#932667FF", #viridisLite::inferno(1, begin = 0.4),
        high = "#F6D645FF", #inferno(1, begin = 0.9),
        midpoint = 0,
        mid = "#000004FF", #inferno(1, begin = 0)
      ) +
      theme(
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = hj),
        axis.title.x = element_text(hjust = hj),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points")
      )
    gcPlot <- .updateThemeFromDots(gcPlot, ...)

    ## Get the PWF status
    status <- getSummary(x)
    status <- subset(status, Category == "Per sequence GC content")
    status$Filename <- factor(
      labels[status$Filename], levels = levels(df$Filename)
    )
    status <- dplyr::right_join(
      status, unique(df["Filename"]), by = "Filename"
    )

    if (missing(pwfCols)) pwfCols <- ngsReports::pwf
    gcPlot <-
      .prepHeatmap(gcPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

  }

  gcPlot
}
)
#' @importFrom dplyr bind_rows
#' @importFrom scales percent
#' @rdname plotGcContent-methods
#' @export
setMethod(
  "plotGcContent", signature = "FastpData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".fastp.*", theoreticalGC = TRUE,
    gcType = c("Genome", "Transcriptome"), species = "Hsapiens", GCobject,
    Fastafile, n = 1e6, plotType = "bar", fillCol = c("navyblue", "red3"), ...
  ){

    plotType <- match.arg(plotType)
    mod <- "Summary"
    steps <-  c("Before_filtering", "After_filtering")
    data <- getModule(x, mod)
    data <- data[names(data) %in% steps]
    df <- bind_rows(data, .id = "Step")
    cols <- c("Step", "Filename", "gc_content")
    df <- df[names(df) %in% cols]
    df$Step <- factor(df$Step, levels = steps)

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, pattern = pattern,...)
    labels <- labels[names(labels) %in% df$Filename]
    df$Filename <- factor(labels[df$Filename], levels = labels)

    ## Get the theoretical mean GC
    exp_mean <- c()
    if (theoreticalGC) {
      if (!missing(Fastafile)) {
        rdLen <- data[[1]]$read1_mean_length
        gc_df <- estGcDistn(Fastafile, n, rdLen)
        gc_df$GC_Content
        species <- "Freq"
      }
      else{
        if (missing(GCobject)) GCobject <- ngsReports::gcTheoretical
        gcType <- stringr::str_to_title(gcType)
        gcType <- match.arg(gcType)
        avail <- gcAvail(GCobject, gcType)
        species <- match.arg(species, avail$Name)
        gc_df <- getGC(GCobject, species, gcType)
      }
      gc_df$p <- cumsum(gc_df[[species]])
      gc_df$d <- abs(gc_df$p - 0.5)
      gc_df <- gc_df[which(rank(gc_df$d) <= 2),]
      exp_mean <- weighted.mean(gc_df$GC_Content / 100, 1 / gc_df$d)
    }

    fillCol <- rep_len(fillCol, length(steps))
    df[["% GC"]] <- percent(df$gc_content, 0.1)
    p <- NULL
    if (plotType == "bar"){

      p <- ggplot(
        df,
        aes(
          !!sym("Step"), !!sym("gc_content"), fill = !!sym("Step"),
          main = Filename, label = !!sym("% GC")
        )
      ) +
        geom_col() +
        geom_hline(yintercept = exp_mean) +
        scale_fill_manual(values = fillCol) +
        scale_y_continuous(labels = percent, expand = expansion(c(0, 0.05))) +
        labs(y = "% GC") +
        ggtitle(unique(df$Filename)) +
        theme_bw()

      if (usePlotly) {
        p <- p + theme(legend.position = "none")
        hv <- c("x", "main", "label")
        p <- plotly::ggplotly(p, tooltip = hv)
      }
    }

    p
  }
)
#' @importFrom dplyr bind_rows
#' @importFrom scales percent
#' @rdname plotGcContent-methods
#' @export
setMethod(
  "plotGcContent", signature = "FastpDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".fastp.*", theoreticalGC = TRUE,
    gcType = c("Genome", "Transcriptome"), species = "Hsapiens", GCobject,
    Fastafile, n=1e6, plotType = "bar", fillCol = c("navyblue", "red3"), ...
  ){

    plotType <- match.arg(plotType)
    mod <- "Summary"
    steps <-  c("Before_filtering", "After_filtering")
    data <- getModule(x, mod)
    data <- data[names(data) %in% steps]
    df <- bind_rows(data, .id = "Step")
    cols <- c("Step", "Filename", "gc_content")
    df <- df[names(df) %in% cols]
    df$Step <- factor(df$Step, levels = steps)

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, pattern = pattern,...)
    labels <- labels[names(labels) %in% df$Filename]
    df$Filename <- factor(labels[df$Filename], levels = labels)

    ## Get the theoretical mean GC
    exp_mean <- c()
    if (theoreticalGC) {
      if (!missing(Fastafile)) {
        rdLen <- data[[1]]$read1_mean_length
        gc_df <- estGcDistn(Fastafile, n, rdLen)
        gc_df$GC_Content
        species <- "Freq"
      }
      else{
        if (missing(GCobject)) GCobject <- ngsReports::gcTheoretical
        gcType <- stringr::str_to_title(gcType)
        gcType <- match.arg(gcType)
        avail <- gcAvail(GCobject, gcType)
        species <- match.arg(species, avail$Name)
        gc_df <- getGC(GCobject, species, gcType)
      }
      gc_df$p <- cumsum(gc_df[[species]])
      gc_df$d <- abs(gc_df$p - 0.5)
      gc_df <- gc_df[which(rank(gc_df$d) <= 2),]
      exp_mean <- weighted.mean(gc_df$GC_Content / 100, 1 / gc_df$d)
    }

    if (plotType == "bar") {

      ## Add columns for nicer plotting with plotly
      df[["% GC"]] <- percent(df$gc_content, 0.01)
      fillCol <- rep_len(fillCol, length(steps))
      names(fillCol) <- steps

      p <- ggplot(
        df,
        aes(
          Filename, !!sym("gc_content"), fill = !!sym("Step"),
          rate = !!sym("% GC")
        )
      ) +
        geom_col(position = "dodge") +
        geom_hline(yintercept = exp_mean) +
        scale_y_continuous(expand = expansion(c(0, 0.05)), labels = percent) +
        scale_fill_manual(values = fillCol) +
        labs(y = "GC Content (%)") +
        theme_bw()

      p <- .updateThemeFromDots(p, ...)
    }

    if (usePlotly) {
      p <- p + theme(legend.position = "none")
      hv <- c("Filename", "% GC", "Step")
      p <- plotly::ggplotly(p, tooltip = hv)
    }

    p
  }
)
