#' @title Plot the PASS/WARN/FAIL information
#'
#' @description Extract the PASS/WARN/FAIL summaries and plot them
#'
#' @details This uses the standard ggplot2 syntax to create a three colour plot.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly \code{logical}. Generate an interactive plot using plotly
#' @param ... Used to pass various potting parameters to theme.
#' @param gridlineWidth,gridlineCol Passed to geom_hline and geom_vline to determine
#' width and colour of gridlines
#'
#' @return A ggplot2 object (\code{usePlotly = FALSE})
#' or an interactive plotly object (\code{usePlotly = TRUE})
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
#' # Check the overall PASS/WARN/FAIL status
#' plotSummary(fdl)
#'
#' # Interactive plot
#' plotSummary(fdl, usePlotly = TRUE)
#'
#' @importFrom dplyr bind_cols
#' @importFrom magrittr %>%
#' @importFrom magrittr set_names
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom grid unit
#'
#' @name plotSummary
#' @rdname plotSummary-methods
#' @export
setGeneric("plotSummary",function(x, usePlotly = FALSE, ...){standardGeneric("plotSummary")})
#' @aliases plotSummary,character
#' @rdname plotSummary-methods
#' @export
setMethod("plotSummary", signature = "character",
          function(x, usePlotly = FALSE, ...){
            if (length(x) == 1) stop("plotSummary() can only be called on two or more files.")
            x <- getFastqcData(x)
            plotSummary(x, usePlotly,...)
          }
)
#' @aliases plotSummary,FastqcFileList
#' @rdname plotSummary-methods
#' @export
setMethod("plotSummary", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSummary(x, usePlotly,...)
          }
)
#' @aliases plotSummary,FastqcDataList
#' @rdname plotSummary-methods
#' @export
setMethod("plotSummary", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, pwfCols, ..., gridlineWidth = 0.2, gridlineCol = "grey20"){

            df <- getSummary(x)

            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))
            fillCol <- getColours(pwfCols)

            df$Category <- factor(df$Category, levels = unique(df$Category))
            df$Status <- factor(df$Status, levels = rev(c("PASS", "WARN", "FAIL")))

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Add any new labels
            df$Filename <- labels[df$Filename]
            df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

            if(usePlotly){
              nx <- length(x)
              ny <- length(unique(df$Category))

              df$StatusNum <- as.integer(df$Status)
              sumPlot <- ggplot(df, aes_string(x = "Category", y = "Filename", fill = "StatusNum", key = "Status")) +
                geom_tile(colour = gridlineCol) +
                geom_vline(xintercept = seq(1.5, nx), colour = gridlineCol, size = gridlineWidth) +
                geom_hline(yintercept = seq(1.5, ny), colour = gridlineCol, size = gridlineWidth) +
                scale_fill_gradientn(colours = c(fillCol["FAIL"],
                                                 fillCol["WARN"],
                                                 fillCol["PASS"]),
                                     values = c(0,1)) +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                # labs(x="QC Category", y="Filename") +
                theme_bw() +
                theme(axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      axis.title = element_blank(),
                      plot.margin = unit(c(0.01, 0.01, 0.01, 0.04), "npc"),
                      legend.position = "none")

              # Add any parameters from dotArgs
              if (!is.null(userTheme)) sumPlot <- sumPlot + userTheme
              suppressMessages(ggplotly(sumPlot, tooltip = c("Category", "Filename", "Status")))
            }
            else{

              sumPlot <- ggplot(df, aes_string(x = "Category", y = "Filename", fill = "Status")) +
                geom_tile(colour = gridlineCol, size = gridlineWidth) +
                scale_fill_manual(values = fillCol) +
                labs(x="QC Category", y="Filename") +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

              # Add any parameters from dotArgs
              if (!is.null(userTheme)) sumPlot <- sumPlot + userTheme
              sumPlot
            }

          }
)
