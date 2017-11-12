#' @title Draw an N Content Plot
#'
#' @description Draw an N Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the N_Content from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' When \code{x} is a single FastqcFile, or FastqcData object line plots will always
#' be drawn for all Ns.
#' Otherwise, users can select line plots or heatmaps.
#'
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical}. Output as ggplot2 (default) or plotly object.
#' @param warn,fail The default values for warn and fail are 5 and 10 respectively (i.e. precentages)
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default
#' @param lineCol Defaults to red
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A standard ggplot2 object, or an interactive plotly object
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
#' # The default plot
#' plotNContent(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fileName(fdl))
#' plotNContent(fdl[r1])
#'
#' # Plot just the Universal N
#' # and change the y-axis using ggplot2::scale_y_continuous
#' library(ggplot2)
#' plotNContent(fdl, NType ="Universal", plotType = "line") +
#' scale_y_continuous()
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_colour_discrete
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom dplyr vars
#' @importFrom dplyr data_frame
#' @importFrom dplyr funs
#' @importFrom magrittr %>%
#'
#' @name plotNContent
#' @rdname plotNContent-methods
#' @export
setGeneric("plotNContent",function(x, usePlotly = FALSE, ...){standardGeneric("plotNContent")})
#' @aliases plotNContent,character
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotNContent(x, usePlotly,...)
          }
)
#' @aliases plotNContent,FastqcFile
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotNContent(x, usePlotly,...)
          }
)
#' @aliases plotNContent,FastqcFileList
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotNContent(x, usePlotly,...)
          }
)
#' @aliases plotNContent,FastqcData
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, warn = 5, fail = 20, pwfCols, ..., lineCol = "red"){

            # Get the NContent
            df <- Per_base_N_content(x)
            colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

            # Return an empty plot if required
            # if (sum(df$Percentage) == 0){
            #   label <- paste0(fileName(x), ": No N content detected")
            #   return(emptyPlot(label))
            # }

            # Sort out the colours
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))
            pwfCols <- setAlpha(pwfCols, 0.2)

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated Totals")

            df$Filename <- labels[df$Filename]
            df$Base <- factor(df$Base, levels = unique(df$Base))
            df$xValue <- as.integer(df$Base)

            # Setup the BG colours
            rects <- dplyr::data_frame(xmin = 0,
                                       xmax = max(df$xValue),
                                       ymin = c(0, warn, fail),
                                       ymax = c(warn, fail, 100),
                                       Status = c("PASS", "WARN", "FAIL"))

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            nPlot <- ggplot(df) +
              geom_rect(data = rects,
                        aes_string(xmin = "xmin", xmax = "xmax",
                                   ymin = "ymin", ymax = "ymax", fill = "Status")) +
              geom_line(aes_string(x = "xValue", y = "Percentage"), colour = lineCol) +
              scale_fill_manual(values = getColours(pwfCols)) +
              scale_x_continuous(breaks = unique(df$xValue), labels = levels(df$Base), expand = c(0,0)) +
              scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
              facet_wrap(~Filename) +
              labs(x = "Position in Read",
                   y = "%N") +
              guides(fill = FALSE) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

            # Add the basic customisations
            if (!is.null(userTheme)) nPlot <- nPlot + userTheme

            if (usePlotly){
              nPlot <- suppressMessages(
                plotly::ggplotly(nPlot + theme(legend.position = "none"),
                                 hoverinfo = c("x", "y", "colour"))
              )
              # Set the hoverinfo for bg rectangles to the vertices only,
              # This will effectively hide them
              nPlot$x$data[[1]]$hoveron <- "points"
              nPlot$x$data[[2]]$hoveron <- "points"
              nPlot$x$data[[3]]$hoveron <- "points"
            }

            nPlot

          }
)
