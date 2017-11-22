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
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
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
#' plotNContent(fdl[[1]])
#'
#'
#' @importFrom dplyr vars funs
#' @import ggplot2
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
              nPlot <- nPlot + xlab("")
              nPlot <- suppressMessages(
                plotly::ggplotly(nPlot + theme(legend.position = "none"),
                                 hoverinfo = c("x", "y", "colour")) %>%
                  plotly::layout(xaxis = list(title = "Position in Read (bp)"), margin = list(b = 70))
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
#' @aliases plotNContent,FastqcDataList
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, warn = 5, fail = 20, pwfCols,
                   clusterNames = FALSE, dendrogram = FALSE, ...){

            # Get the NContent
            df <- Per_base_N_content(x)
            colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

            # Sort out the colours
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated Totals")

            ## fill bins up to the max sequence length
            df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))
            df <- split(df, f = df$Filename) %>%
              lapply(function(x){
                Longest_sequence <- max(x$Start)
                dfFill <- data.frame(Start = 1:Longest_sequence)
                x <- dplyr::right_join(x, dfFill, by = "Start") %>%
                  zoo::na.locf()
              }) %>%
              dplyr::bind_rows()

            # get ready for clustering and making the key
            # this is the simplest way, in my mind, to do the key
            df$Start <- as.integer(df$Start)
            df <- df[colnames(df) %in% c("Filename", "Start", "Percentage")]
            df <- reshape2::dcast(df, Filename ~ Start, value.var = "Percentage")

            #cluster
            if(clusterNames){
              xx <- df[!colnames(df) == "Filename"]
              xx[is.na(xx)] <- 0
              clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
              row.ord <- order.dendrogram(clus)
              df <- df[row.ord,]
            }

            key <- df$Filename
            df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Start", value.name = "Percentage")
            df$Filename <- labels[df$Filename]

            # Reverse the factor levels for a better looking default plot
            df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))
            df$Percentage <- as.numeric(df$Percentage)
            df$Start <- as.integer(as.character(df$Start))


            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Reverse the filename levels for an alphabetic plot
            df$Filename <- factor(df$Filename, levels = unique(df$Filename))

            # Return an empty plot if required
            allZero <- ifelse(sum(df$Percentage, na.rm = TRUE) == 0, TRUE, FALSE)

            if (allZero){
              # will put the message only in the center of the plot
              label_df <- dplyr::data_frame(
                Filename = levels(df$Filename)[floor(mean(as.integer(df$Filename)))],
                Start = df$Start,
                text = "No N Content Detected")
              nPlot <- ggplot(df) +
                geom_blank(aes_string("Start", "Filename")) +
                geom_text(data = label_df, aes_string(label = "text"), x = max(df$Start)/2, y = length(x)/2) +
                scale_x_continuous(expand = c(0,0)) +
                scale_y_discrete(expand = c(0, 0)) +
                labs(x = "Position in Read",
                     y = "Filename",
                     fill = "%N") +
                theme_bw() +
                theme(panel.background = element_rect(fill = "white"),
                      panel.grid = element_blank())
            }
            else{
              nPlot <- ggplot(df, aes_string("Start", "Filename", fill = "Percentage")) +
                geom_tile() +
                scale_fill_pwf(df$Percentage, pwfCols, breaks = c(0, warn, fail, 101), passLow = TRUE, na.value = "white") +
                scale_x_continuous(expand = c(0,0)) +
                scale_y_discrete(expand = c(0, 0)) +
                labs(x = "Position in Read",
                     y = "Filename",
                     fill = "%N") +
                theme_bw() +
                theme(panel.background = element_blank())
            }

            # Add the basic customisations
            if (!is.null(userTheme)) nPlot <- nPlot + userTheme

            if (usePlotly){
              # Reset the status using current values
              status <- dplyr::summarise_at(dplyr::group_by(df, Filename),
                                            vars("Percentage"), funs(Percentage = max), na.rm = TRUE)
              status$Status <- cut(status$Percentage, breaks = c(0, warn, fail, 101), include.lowest = TRUE,
                                   labels = c("PASS", "WARN", "FAIL"))

              # Form the sideBar for each adapter
              sideBar <- makeSidebar(status, key, pwfCols = pwfCols)

              # Customise for plotly
              nPlot <- nPlot +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
              if (!is.null(userTheme)) nPlot <- nPlot + userTheme

              if (clusterNames && dendrogram){
                dx <- ggdendro::dendro_data(clus)
                dendro <- ggdend(dx$segments) +
                  coord_flip() +
                  scale_y_reverse(expand = c(0, 0)) +
                  scale_x_continuous(expand = c(0, 0.5))

                nPlot <- suppressWarnings(
                  suppressMessages(
                    plotly::subplot(dendro, sideBar, nPlot, widths = c(0.1,0.08,0.82),
                                    margin = 0.001, shareY = TRUE)
                  ))
              }
              else{

                nPlot <- suppressWarnings(suppressMessages(
                  plotly::subplot(plotly::plotly_empty(), sideBar, nPlot, widths = c(0.1,0.08,0.82), margin = 0.001, shareY = TRUE) %>%
                    plotly::layout(annotations = list(text = "Filename", showarrow = FALSE,
                                                      textangle = -90))
                ))
              }
              nPlot <- nPlot %>%
                plotly::layout(xaxis3 = list(title = "Position in Read (bp)"), margin = list(b = 45))

            }

            nPlot

          }
)
