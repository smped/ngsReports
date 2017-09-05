#' @title Plot the combined Sequence_Duplication_Levels information
#'
#' @description Plot the Sequence_Duplication_Levels information for a set of FASTQC reports
#'
#' @details
#' This extracts the Sequence_Duplication_Levels from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param warn,fail The default values for warn and fail are 20 and 50 respectively (i.e. percentages)
#' @param lineCols Colours of the lines drawn for individual libraries
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#'
#' @return A standard ggplot2 or plotly object
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
#' # Draw the default plot for a single file
#' plotDuplicationLevels(fdl[[1]])
#'
#' @name plotDuplicationLevels
#' @rdname plotDuplicationLevels-methods
#' @export
setGeneric("plotDuplicationLevels",function(x, usePlotly = FALSE, ...){standardGeneric("plotDuplicationLevels")})
#' @aliases plotDuplicationLevels,character
#' @rdname plotDuplicationLevels-methods
#' @export
setMethod("plotDuplicationLevels", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotDuplicationLevels(x, usePlotly,...)
          }
)
#' @aliases plotDuplicationLevels,FastqcFile
#' @rdname plotDuplicationLevels-methods
#' @export
setMethod("plotDuplicationLevels", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotDuplicationLevels(x, usePlotly,...)
          }
)
#' @aliases plotDuplicationLevels,FastqcData
#' @rdname plotDuplicationLevels-methods
#' @export
setMethod("plotDuplicationLevels", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 20, fail = 50, lineCols = c("red", "blue"), ...){

            df <- Sequence_Duplication_Levels(x)
            df <- reshape2::melt(df, id.vars = c("Filename", "Duplication_Level"),
                                 value.name = "Percentage", variable.name = "Type")
            df$Duplication_Level <- factor(df$Duplication_Level, levels = unique(df$Duplication_Level))
            df$Type <- gsub("Percentage_of_", "% ", df$Type)
            df$Type <- stringr::str_to_title(df$Type)
            df$Type <- paste(df$Type, "sequences")
            df

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
            df$Filename <- labels[df$Filename]
            df$x <- as.integer(df$Duplication_Level)
            df$Percentage <- round(df$Percentage, 2)

            # Get any theme arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Sort out the colours
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))
            pwfCols <- setAlpha(pwfCols, 0.2)
            rects <- dplyr::data_frame(xmin = min(df$x) - 1,
                                       xmax = max(df$x) + 1,
                                       ymin = c(0, 100 - fail, 100 - warn),
                                       ymax = c(100 - fail, 100 - warn, 100),
                                       Status = c("FAIL", "WARN", "PASS"))

            dupPlot <- ggplot(data = df) +
              geom_rect(data = rects, aes_string(xmin = "xmin", xmax = "xmax",
                                                 ymin = "ymin", ymax = "ymax",
                                                 fill = "Status")) +
              geom_line(aes_string(x = "x", y = "Percentage", colour = "Type", group = "Type")) +
              scale_fill_manual(values = getColours(pwfCols)) +
              scale_colour_manual(values = lineCols) +
              scale_x_continuous(breaks = unique(df$x), labels = levels(df$Duplication_Level), expand = c(0, 0)) +
              scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
              labs(x = "Sequence Duplication Level",
                   colour = c()) +
              guides(fill = FALSE) +
              ggtitle(unique(df$Filename)) +
              theme_bw() +
              theme(legend.position = c(1, 1),
                    legend.justification = c(1, 1),
                    legend.background = element_rect(colour = "black", size = 0.2),
                    plot.title = element_text(hjust = 0.5))
            if (!is.null(userTheme)) dupPlot <- dupPlot + userTheme

            # browser()

            if (usePlotly){

              dupPlot <- dupPlot +
                xlab("") +
                theme(legend.position = "none")
              dupPlot <- suppressMessages(
                plotly::ggplotly(dupPlot, tooltip = c("colour", "Percentage"))
              )
              dupPlot$x$data[[1]]$hoveron <- "points"
              dupPlot$x$data[[2]]$hoveron <- "points"
              dupPlot$x$data[[3]]$hoveron <- "points"
            }

            dupPlot

          }
)
# plotDuplicationLevels <- function(x, subset, type = ".+",
#                                   trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",
#                                   usePlotly = FALSE){
#
#   stopifnot(grepl("(Fastqc|character)", class(x)))
#
#   if (missing(subset)){
#     subset <- rep(TRUE, length(x))
#   }
#   stopifnot(is.logical(subset))
#   stopifnot(length(subset) == length(x))
#   stopifnot(is.logical(trimNames))
#   stopifnot(is.character(type))
#
#   x <- x[subset]
#   df <- tryCatch(Sequence_Duplication_Levels(x))
#
#   # Check the pattern contains a capture
#   if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
#     df$Filename <- gsub(pattern[1], "\\1", df$Filename)
#     # These need to be checked to ensure non-duplicated names
#     if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
#   }
#
#     # Melt for a convenient long form
#   df <- reshape2::melt(df, id.vars = c("Filename", "Duplication_Level"),
#                        variable.name = "Type", value.name = "Percent")
#
#   # Ensure the Duplication_level plots in the correct order
#   df$Duplication_Level <- factor(df$Duplication_Level,
#                                  levels = c(1:9, paste0(">", c(10, 50, 100, 500, "1k", "5k", "10k+"))))
#
#   # Tidy the output to look like the original report
#   df$Type <- gsub("Percentage_of_", "% ", df$Type)
#   df$Type <- stringr::str_to_title(df$Type)
#   df$Type <- paste(df$Type, "sequences")
#
#   # Restrict to a given type if requested
#   df <- dplyr::filter(df, grepl(type, Type))
#   if(nrow(df) == 0) stop("Invalid type selected")
#
#   dupPlot <- ggplot(df, aes(x = as.integer(Duplication_Level), y = Percent)) +
#     annotate("rect",
#                       xmin = seq(1.5, 17, by = 2), xmax = seq(2.5, 17, by = 2),
#                       ymin = 0, ymax = Inf, fill = "grey80", alpha = 0.5) +
#     geom_line(aes(colour = Filename, group = Filename)) +
#     facet_wrap(~Type, ncol = 1) +
#     scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
#                                 breaks = seq(0, 100, by = 20)) +
#     scale_x_continuous(breaks = seq_along(levels(df$Duplication_Level)),
#                                 labels = levels(df$Duplication_Level),
#                                 limits = c(0, length(levels(df$Duplication_Level))) + 0.5,
#                                 expand = c(0, 0)) +
#     labs(x = "Duplication Level", y = "Percent (%)") +
#     theme_bw() +
#     theme(panel.grid.major.x = element_blank(),
#                    panel.grid.minor.x = element_blank())
#
#     if(usePlotly){
#
#       if(length(x) == 1) dupPlot <- dupPlot + labs(colour = "")
#
#       dupPlot <- ggplotly(dupPlot)  %>% layout(legend = list(orientation = "h"))
#
#     }
#
#     # And draw the plot
#     dupPlot
# }
