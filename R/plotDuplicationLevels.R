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
#' @param deduplication Plot Duplication levels 'pre' or 'post' deduplication. Can only take values "pre" and "post"
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
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
#' plotDuplicationLevels(fdl)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 ggtitle
#'
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
#' @aliases plotDuplicationLevels,FastqcFileList
#' @rdname plotDuplicationLevels-methods
#' @export
setMethod("plotDuplicationLevels", signature = "FastqcFileList",
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
#' @aliases plotDuplicationLevels,FastqcDataList
#' @rdname plotDuplicationLevels-methods
#' @export
setMethod("plotDuplicationLevels", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, deduplication = "pre", clusterNames = FALSE, dendrogram = FALSE, labels, pwfCols, ...){

            df <- tryCatch(Sequence_Duplication_Levels(x))
            stopifnot(deduplication %in% c("pre", "post"))
            type <- c(pre = "Percentage_of_total", post = "Percentage_of_deduplicated")[deduplication]
            df <- df[c("Filename", "Duplication_Level",type)]
            df[[type]] <- round(df[[type]], 2)

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
            df$Filename <- labels[df$Filename]

            # Get any theme arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            if(clusterNames){
              # Grab the main columns & cast from long to wide
              mat <- reshape2::acast(df[c("Filename", "Duplication_Level", type)], Filename ~ Duplication_Level, value.var = type)
              mat[is.na(mat)] <- 0
              clus <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
              row.ord <- order.dendrogram(clus)
              mat <- mat[row.ord,]
              o <- rownames(mat)
            }
            else{
              o <- unique(df$Filename)
            }
            df$Filename <- factor(df$Filename, levels = o)
            df$Duplication_Level <- factor(df$Duplication_Level, levels = unique(df$Duplication_Level))
            df$x <- as.integer(df$Duplication_Level)
            key <- unique(df$Filename)

            fillLabel <- gsub("Percentage_of_", "% ", type)
            fillLabel <- paste(fillLabel, "\nSequences")
            fillLabel <- stringr::str_to_title(fillLabel)

            dupPlot <- ggplot(df, aes_string("x", "Filename")) +
              geom_tile(aes_string(fill = type)) +
              scale_x_continuous(breaks = unique(df$x), labels = levels(df$Duplication_Level), expand = c(0, 0)) +
              scale_y_discrete(expand = c(0, 0)) +
              scale_fill_gradientn(colours = inferno(50)) +
              labs(x = "Sequence Duplication Level",
                   fill = fillLabel) +
              theme(panel.grid = element_blank(),
                    panel.background = element_blank())

            if (usePlotly){

            dupPlot <- dupPlot +
              theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())

              t <- getSummary(x)
              t <- t[t$Category == "Sequence Duplication Levels",]
              t$Filename <- factor(labels[t$Filename], levels = levels(df$Filename))
              t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")

              if (missing(pwfCols)) pwfCols <- ngsReports::pwf

              sideBar <- ggplot(t, aes(x = 1, y = Filename, key = key)) +
                geom_tile(aes_string(fill = "Status")) +
                scale_fill_manual(values = getColours(pwfCols)) +
                theme(panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      legend.position="none",
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank())
              sideBar <- suppressMessages(plotly::ggplotly(sideBar, tooltip = c("Status", "Filename")))

              #plot dendrogram
              if(dendrogram && clusterNames){

                dx <- ggdendro::dendro_data(clus)
                dendro <- ggdend(dx$segments) +
                  coord_flip() +
                  scale_y_reverse(expand = c(0, 1)) +
                  scale_x_continuous(expand = c(0,1))

                dupPlot <- suppressWarnings(
                  suppressMessages(
                    plotly::subplot(dendro, sideBar, dupPlot, widths = c(0.1,0.1,0.8),
                                    margin = 0, shareY = TRUE) %>%
                      plotly::layout(xaxis3 = list(title = "Sequence Duplication Levels"))
                  ))
              }
              else{
                dupPlot <- suppressWarnings(
                  suppressMessages(
                    plotly::subplot(plotly::plotly_empty(), sideBar, dupPlot,
                                    widths = c(0.1,0.1,0.8), margin = 0, shareY = TRUE) %>%
                      plotly::layout(xaxis3 = list(title = "Sequence Duplication Levels"),
                                     annotations = list(text = "Filename", showarrow = FALSE,
                                                        textangle = -90))
                  )
                )
              }

            }

            dupPlot

          }
)
