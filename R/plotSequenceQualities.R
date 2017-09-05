#' @title Plot the Per Sequence Quality Scores
#'
#' @description Plot the Per Sequence Quality Scores for a set of FASTQC reports
#'
#' @details Plots the distribution of average sequence quality scores across the set of files.
#' Values can be plotted either as counts (\code{counts = TRUE}) or as frequencies (\code{counts = FALSE}).
#'
#' Any faceting or scale adjustment can be performed after generation of the initial plot,
#' using the standard methods of ggplot2 as desired.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param warn,fail The default values for WARN and FAIL in FASTQC
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param ... Used to pass additional attributes to theme()
#' @param expand the default expansion of the x-axis. Passed to scale_x_continuous()
#' @param alpha Passed to background colours
#'
#' @return A ggplot2 object or interactive plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # Draw the defualt plot
#' plotSequenceQualities(fdl)
#'
#' # Customise using the R1 subset, plotting counts,
#' # and faceting after the initial function call
#' library(ggplot2)
#' plotSequenceQualities(fdl, counts = TRUE) +
#'   facet_wrap(~Filename, ncol = 2, scales = "free_y")
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @importFrom plotly ggplotly
#' @importFrom plotly add_trace
#' @importFrom grDevices adjustcolor
#' @importFrom magrittr %>%
#'
#' @export
plotSequenceQualities <- function(x, usePlotly = FALSE, labels, counts = FALSE, warn = 30, fail = 20, pwfCols,
                                  ..., expand = c(0.02, 0), alpha = 0.2){

  df <- tryCatch(Per_sequence_quality_scores(x))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))
  cols <- getColours(pwfCols)

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  minQ <- min(df$Quality)

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  if (!counts){
    Count <- NULL # To avoid NOTE messages in R CMD check
    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>% dplyr::mutate(Freq = Count / sum(Count))
    df <- dplyr::ungroup(df)
    qualPlot <- ggplot(df, aes_string(x = "Quality", y = "Freq", colour = "Filename")) +
      ylab("Frequency")

  }
  else{

    # Initialise the plot using counts
    qualPlot <- ggplot(df, aes_string(x = "Quality", y = "Count", colour = "Filename")) +
      ylab("Count")

  }

  qualPlot <- qualPlot +
    annotate("rect", xmin = warn, xmax = Inf, ymin = -Inf, ymax = Inf, fill = cols["PASS"],
             alpha = alpha) +
    annotate("rect", xmin = fail, xmax = warn, ymin = -Inf, ymax = Inf, fill = cols["WARN"],
             alpha = alpha) +
    annotate("rect", xmin = -Inf, xmax = fail, ymin = -Inf, ymax = Inf, fill = cols["FAIL"],
             alpha = alpha) +
    scale_x_continuous(limits = c(minQ, 40), breaks = seq(0, 40, by = 10), expand = expand) +
    scale_colour_discrete(labels = labels) +
    geom_line() +
    xlab("Mean Sequence Quality Per Read (Phred Score)") +
    theme_bw()

  if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

  if(usePlotly){
    cutOffs <- data.frame(pass = warn, Filename = df$Filename, warn = fail, fail = 0, top =  41)
    ymax <- 1
    if (counts) ymax <- max(df$Count)

    qualPlot <- qualPlot + labs(colour = "")
    qualPlot <- suppressMessages(ggplotly(qualPlot))
    # These are still highly problematic.
    # I wonder why layout using shape doesn't work?
    qualPlot <- suppressWarnings(
      add_trace(qualPlot,
                data = cutOffs, x = ~top, type = 'scatter', mode = 'lines',
                line = list(color = NULL),
                showlegend = FALSE, fillopacity = alpha,
                xmin = 0, xmax = Inf, ymin = 0, ymax =  Inf,
                hoverinfo = "none"))
    qualPlot <- suppressWarnings(
      add_trace(qualPlot,
                data = cutOffs, x = ~pass, type = 'scatter', mode = 'lines',
                fill = 'tonexty', fillcolor=adjustcolor(cols["PASS"], alpha.f = alpha),
                line = list(color = adjustcolor(cols["PASS"], alpha.f = alpha)),
                showlegend = FALSE,
                xmin = 0, xmax = Inf, name = "PASS", ymin = warn, ymax = 0,
                hoverinfo = "none"))
    qualPlot <- suppressWarnings(
      add_trace(qualPlot,
                data = cutOffs, x = ~warn, type = 'scatter', mode = 'lines',
                fill = 'tonexty', fillcolor=adjustcolor(cols["WARN"], alpha.f = alpha),
                line = list(color = adjustcolor(cols["WARN"], alpha.f = alpha)),
                showlegend = FALSE,
                xmin = 0, xmax = Inf, name = "WARN", ymin = -Inf, ymax = 0,
                hoverinfo = "none"))
    qualPlot <- suppressWarnings(
      add_trace(qualPlot,
                data = cutOffs, x = ~fail, type = 'scatter', mode = 'lines',
                fill = 'tonexty', fillcolor=adjustcolor(cols["FAIL"], alpha.f = alpha),
                line = list(color = adjustcolor(cols["FAIL"], alpha.f = alpha)),
                showlegend = FALSE,
                xmin = -Inf, xmax = Inf,
                name = "FAIL", ymin = -Inf, ymax = 0, hoverinfo = "none"))

  }

  # Draw the plot
  qualPlot

}
