#' @title Plot Overrepresented Kmers
#'
#' @description Plot Overrepresented Kmers
#'
#' @details As the Kmer Content module present in FastQC reports is relatively
#' uninformative, and omitted by default in later versions of FastQC, these
#' are rudimentary plots.
#'
#' Plots for `FastqcData` objects replicate those contained in a FastQC
#' report, whilst the heatmap generated from `FastqcDataList` objects
#' simply show the location and abundance of over-represented Kmers.
#'
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param n `numeric`. The number of Kmers to show.
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' @param pattern regex to drop from the end of filenames
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#' @param linewidth Passed to `geom_line()`
#' @param pal The colour palette. If the vector supplied is less than n,
#' `grDevices::colorRampPalette()` will be used
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail values in the plot
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param heatCol Colour palette used for the heatmap. Default is `inferno`
#' from the viridis set of palettes
#' @param heat_w Relative width of any heatmap plot components
#' @param module The module to obtain data from when using a FastpData object
#' @param reads Either read1 or read2. Only used when using a FastpData object
#' @param fill A ggplot2 continuous scale to be used for the colour palette
#' @param trans Function for transforming the count/mean ratio. Set as NULL
#' to use the ratio without transformation
#' @param labelSize,labelCol Size & colour for the kmer labels drawn in
#' each each cell
#' @param readsBy Strategy for visualising both read1 and read2. Can be set to
#' show each set of reads by facet, or within the same plot taking the mean of
#' the enrichment above mean, or the difference in the enrichment above mean
#' @param fill_lab Label for fill legend. Set to NULL to have no label.
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
#' plotKmers(fdl[[1]])
#'
#' # Use a FastpData object
#' fl <- system.file("extdata", "fastp.json", package = "ngsReports")
#' fp <- FastpData(fl)
#' plotKmers(fp)
#'
#'
#' @docType methods
#'
#' @importFrom dplyr desc
#' @importFrom grDevices hcl.colors
#' @import ggplot2
#' @import tibble
#'
#' @name plotKmers
#' @rdname plotKmers-methods
#' @export
setGeneric("plotKmers", function(x, usePlotly = FALSE, labels, ...){
  standardGeneric("plotKmers")
})
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "ANY", function(
    x, usePlotly = FALSE, labels, ...){
  .errNotImp(x)
}
)
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "character", function(
    x, usePlotly = FALSE, labels, ...){
  x <- FastqcDataList(x)
  if (length(x) == 1) x <- x[[1]]
  plotKmers(x, usePlotly, labels, ...)
}
)
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, pattern = ".fastq.*|.fastp.*|.fq*",
    n = 6, ..., linewidth = 0.5,
    pal = c("red", "blue", "green", "black", "magenta", "yellow")){

  ## Get the basic data frame
  df <- getModule(x, "Kmer_Content")

  if (!length(df)) {
    kMerPlot <- .emptyPlot("No Kmer_Content Module Detected")
    if (usePlotly) kMerPlot <- ggplotly(kMerPlot, tooltip = "")
    return(kMerPlot)
  }

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, pattern = pattern, ...)
  labels <- labels[names(labels) %in% df$Filename]
  df$Filename <- labels[df$Filename]

  ## Get the top kMers
  o <- order(df$PValue, 1/df$Count)
  allK <- unique(df$Sequence)
  n <- tryCatch(as.integer(n))
  n <- min(length(allK), n)
  topK <- unique(df$Sequence[o])[seq_len(n)]

  ## Tidy up the data
  df <- dplyr::filter(df, Sequence %in% topK)
  colnames(df) <- gsub("Max_Obs/Exp_Position", "Base", colnames(df))
  colnames(df) <- gsub("Obs/Exp_Max", "Value", colnames(df))
  df$Position <- gsub("([0-9]*)-[0-9]*", "\\1", df$Base)
  df$Position <- as.integer(df$Position)
  df <- df[c("Filename", "Sequence", "Value", "Position")]

  ## Set the x-axis for plotting
  ## As there will be significant gaps in this data,
  ## the bins used need to be obtained from another slot.
  ## The most complete will be Per_base_sequence_quality
  ## These values can then be incorporated in the final df
  ## for accurate plotting & labelling
  refForX <- unique(getModule(x, "Per_base_sequence_quality")$Base)
  refForX <- tibble(
    Base = as.character(refForX),
    Position = gsub("([0-9]*)-[0-9]*", "\\1", Base)
  )
  refForX$Position <- as.integer(refForX$Position)

  ## In order to get a line plot, zero values need to be added to the
  ## missing positions. The above reference scale for X will be used
  ## to label the missing values
  ## Include all files to ensure all appear in the final plot
  zeros <- expand.grid(
    list(
      Filename = df$Filename,
      Sequence = topK,
      Value = 0,
      Position = refForX$Position
    ),
    stringsAsFactors = FALSE
  )

  ## After the bind_rows, duplicate values will exist at some
  ## positions. Spuriously introduced zeros need to be removed
  df <- dplyr::bind_rows(df, zeros)
  df <- dplyr::arrange(df, Filename, Sequence, Position, desc(Value))
  df <- dplyr::distinct(df, Filename, Sequence, Position, .keep_all = TRUE)
  df <- dplyr::left_join(df, refForX, by = "Position")

  ## Set the Sequence as a factor based on the first position it
  ## appears. This way the colours will appear in order in the guide
  ## as well as the plot
  getKLevels <- function(x){
    x <- subset(x, Value > 0) # Non zero counts
    ## Arrange in order of position
    x <- dplyr::arrange(x, Position, Sequence)
    unique(x$Sequence)
  }
  kMerLevels <- getKLevels(df)
  df$Sequence <- factor(df$Sequence, levels = kMerLevels)
  df <- droplevels(df)

  ## Set the plotting params
  yMax <- max(df$Value)*1.05
  ## Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  ## And the colours
  if (n < length(pal)) {
    pal <- pal[seq_len(n)]
  }
  else {
    pal <- grDevices::colorRampPalette(pal)(n)
  }

  ## Now draw the basic plots
  xLab <- "Position in read (bp)"
  yLab <- expression(paste(log[2], " Obs/Exp"))
  kMerPlot <- ggplot(df, aes(Position, Value, colour = Sequence)) +
    geom_line(linewidth = linewidth) +
    facet_wrap(~Filename) +
    scale_x_continuous(
      breaks = refForX$Position, labels = refForX$Base, expand = c(0.02, 0)
    ) +
    scale_y_continuous(limits = c(0, yMax), expand = c(0, 0)) +
    scale_colour_manual(values = pal) +
    theme_bw() +
    labs(x = xLab, y = yLab,  colour = c())

  ## Check for binned x-axis values to decied whether to rotate
  ## x-axis labels. This should be clear if there are more than
  ## 2 characters in the plotted labels
  binned <- any(grepl("-", df$Base))
  if (binned) {
    kMerPlot <- kMerPlot +
      theme(axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 1
      ))
  }
  if (!is.null(userTheme)) kMerPlot <- kMerPlot + userTheme

  if (usePlotly) {
    yLab <- "Log2 Obs/Exp"
    kMerPlot <- kMerPlot + xlab("")
    kMerPlot <- suppressMessages(
      suppressWarnings(
        plotly::subplot(
          plotly::plotly_empty(),
          kMerPlot,
          widths = c(0.14,0.86)
        )
      )
    )
    kMerPlot <- plotly::layout(
      kMerPlot,
      yaxis2 = list(title = yLab),
      legend = list(orientation = "v", x = 1, y = 1)
    )
  }
  kMerPlot
}
)
#' @rdname plotKmers-methods
#' @export
setMethod(
  "plotKmers", signature = "FastqcDataList",
  function(
    x, usePlotly = FALSE, labels, pattern = ".fastq.*|.fastp.*|.fq*",
    cluster = FALSE, dendrogram = FALSE, pwfCols,
    heatCol = hcl.colors(50, "inferno"), heat_w = 8, ...
  ){

  mod <- "Kmer_Content"
  df <- getModule(x, mod)

  if (!length(df)) {
    kMerPlot <- .emptyPlot("No Overrepresented Kmers Detected")
    if (usePlotly) kMerPlot <- ggplotly(kMerPlot, tooltip = "")
    return(kMerPlot)
  }

  ## Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(.isValidPwf(pwfCols))

  ## Drop the suffix, or check the alternate labels
  labels <- .makeLabels(x, labels, pattern = pattern, ...)
  labels <- labels[names(labels) %in% df$Filename]

  colnames(df) <- gsub("Max_Obs/Exp_Position", "Base", colnames(df))
  colnames(df) <- gsub("Obs/Exp_Max", "Total", colnames(df))
  df$Position <- gsub("([0-9]*)-[0-9]*", "\\1", df$Base)
  df$Position <- as.integer(df$Position)

  ## Now simply add the kmers at each position. Not superinformative.
  df <- dplyr::group_by(df, Filename, Position)
  df <- dplyr::summarise(df, Total = sum(Count))
  df <- dplyr::ungroup(df)

  ## Setup the x-axis
  refForX <- unique(getModule(x, "Per_base_sequence_quality")$Base)
  refForX <- tibble(
    Base = as.character(refForX),
    Position = gsub("([0-9]*)-[0-9]*", "\\1", Base)
  )
  refForX$Position <- as.integer(refForX$Position)

  ## Setup a grid of NAs for positions where no Kmer is detected
  nas <- expand.grid(
    list(
      Filename = unique(df$Filename),
      Total = NA,
      Position = unique(refForX$Position)
    ),
    stringsAsFactors = FALSE
  )
  ## Now remove any spurious NA values
  df <- dplyr::bind_rows(df, nas)
  df <- dplyr::arrange(df, Filename, Position, desc(Total))
  df <- dplyr::distinct(df, Filename, Position, .keep_all = TRUE)

  ## Now define the order for a dendrogram if required
  key <- names(labels)
  clusterDend <- .makeDendro(df, "Filename", "Position", "Total")
  dx <- ggdendro::dendro_data(clusterDend)
  if (dendrogram | cluster) key <- labels(clusterDend)
  if (!dendrogram) dx$segments <- dx$segments[0,]
  ## Now set everything as factors
  df$Filename <- factor(labels[df$Filename], levels = labels[key])

  ## join x axis ticks and df
  df <- dplyr::left_join(df, refForX, by = "Position")

  ## Set up for geom_tile
  df$End <- gsub("[0-9]*-([0-9]*)", "\\1", df$Base)
  df$End <- as.integer(df$End)
  df$Base <- (df$Position + df$End) / 2
  df$Width <- df$End - df$Position
  df$Width[df$Width == 0] <- 1

  xLab <- "Position in Read (bp)"
  hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
  kMerPlot <- ggplot(
    df,
    aes(Base, Filename, fill = Total, width = Width)
  ) +
    geom_tile() +
    ggtitle(gsub("_", " ", mod)) +
    labs(x = xLab, y = c()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), position = "right") +
    scale_fill_gradientn(colors = heatCol, na.value = "white") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"),
      plot.title = element_text(hjust = hj),
      axis.title.x = element_text(hjust = hj)
    )

  dotArgs <- list(...)
  allowed <- names(formals(theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
  if (!is.null(userTheme)) kMerPlot <- kMerPlot + userTheme

  ## Get the PWF status
  status <- getSummary(x)
  status <- status[status$Category == "Kmer Content",]
  status <- subset(status, Filename %in% key)
  status$Filename <- factor(labels[status$Filename], levels = labels[key])

  .prepHeatmap(kMerPlot, status, dx$segments, usePlotly, heat_w, pwfCols)

}
)
#' @importFrom rlang "!!" sym
#' @importFrom forcats fct_rev
#' @importFrom dplyr reframe group_by bind_rows
#' @rdname plotKmers-methods
#' @export
setMethod(
  "plotKmers", signature = "FastpData",
  function(
    x, usePlotly = FALSE, labels, pattern = ".fastq.*|.fastp.*|.fq*",
    module = c("Before_filtering", "After_filtering"),
    reads = c("read1", "read2"), readsBy = c("facet", "mean", "diff"),
    labelSize = 3, labelCol = "grey70", trans = "log2",
    fill = scale_fill_gradient2(),  fill_lab = "Count/\nMean",
    ...
  ){

    ## Check args
    mod <- match.arg(module)
    reads <- match.arg(reads, several.ok = TRUE)
    readsBy <- match.arg(readsBy)
    stopifnot(is(fill, "ScaleContinuous"))

    ## Get the basic data frame
    df <- getModule(x, mod)[reads]
    df <- bind_rows(df, .id = "reads")
    fpName <- gsub(pattern, "", df$Filename)[[1]]
    df[["Filename"]] <- df[["fqName"]]

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, pattern = pattern, ...)
    labels <- labels[names(labels) %in% df$Filename]
    df$Filename <- labels[df$Filename]

    ## Just grab the data we need
    df <- df[c("Filename", "kmer_count")]
    df <- tidyr::unnest(df, !!sym("kmer_count"))
    df[["prefix"]] <- fct_rev(df[["prefix"]])

    ## Make a blank plot if no data is found
    if (!length(df)) {
      msg <- "No kmer Counts Detected"
      p <- .emptyPlot(msg)
      if (usePlotly) p <- ggplotly(p, tooltip = "")
      return(p)
    }

    if (!is.null(trans)) {
      f <- match.fun(trans)
      stopifnot(length(f(1:2)) == 2)
      df[["times_mean"]] <- f(df[["times_mean"]])
      if (!is.null(fill_lab))
        fill_lab <- paste(trans, fill_lab, sep = "\n")
    }

    sub <- c()
    if (readsBy %in% c("mean", "diff")) {
      merge_fun <- match.fun(readsBy)
      cols <- c("kmer", "prefix", "suffix")
      df <- group_by(df, !!!syms(cols))
      df <- reframe(
        df,
        count = merge_fun(!!sym("count")),
        times_mean = merge_fun(!!sym("times_mean"))
      )
      if (!nrow(df)) stop("Please choose another approach for 'readsBy'")
      sub <- paste(
        "Reads combined using",
        ifelse(readsBy == "mean", "mean values", "differences")
      )
    }

    ## Tidy up for plotting
    df[["times_mean"]] <- round(df[["times_mean"]], 3)
    df[["count"]] <- scales::comma(df[["count"]], 1)
    main <- paste0(fpName, ": ", gsub("_", " ", mod))
    p <- ggplot(
      df,
      aes(
        x = !!sym("suffix"), y = !!sym("prefix"), fill = !!sym("times_mean"),
        count = !!sym("count")
      )
    ) +
      geom_raster() +
      geom_text(
        aes(label = !!sym("kmer")), size = labelSize, colour = labelCol
      ) +
      scale_x_discrete(expand = expansion(0), position = "top") +
      scale_y_discrete(expand = expansion(0), position = "left") +
      labs(x = c(), y = c(), fill = fill_lab) +
      ggtitle(main, subtitle = sub) +
      theme(axis.ticks = element_blank()) +
      fill
    if (readsBy == "facet") p <- p + facet_wrap(~Filename, nrow = 1)
    p <- .updateThemeFromDots(p, ...)

    if (usePlotly) {
      p <- p + theme(legend.position = "none")
      hv <- c("kmer", "count", "times_mean")
      p <- suppressWarnings(plotly::ggplotly(p, tooltip = hv))
    }
    p
  }
)
