#' @title Draw an N Content Plot
#'
#' @description Draw an N Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the N_Content from the supplied object and generates a ggplot2
#' object, with a set of minimal defaults.
#' The output of this function can be further modified using the standard 
#' ggplot2 methods.
#'
#' When \code{x} is a single FastqcFile, or FastqcData object line plots will 
#' always be drawn for all Ns.
#' Otherwise, users can select line plots or heatmaps.
#'
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical}. Output as ggplot2 (default) or plotly object.
#' @param warn,fail The default values for warn and fail are 5 and 10 
#' respectively (i.e. percentages)
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours 
#' for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default
#' @param lineCol Defaults to red
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} 
#' then the dendrogram will be displayed.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot
#' plotNContent(fdl[[1]])
#'
#'
#' @import ggplot2
#' @importFrom dplyr vars funs
#' @importFrom zoo na.locf
#'
#' @name plotNContent
#' @rdname plotNContent-methods
#' @export
setGeneric("plotNContent",
           function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
               standardGeneric("plotNContent")
               })
#' @aliases plotNContent,character
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "character",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
              x <- getFastqcData(x)
              plotNContent(x, usePlotly, labels, pwfCols, warn, fail, ...)
          }
)
#' @aliases plotNContent,FastqcFile
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
              x <- getFastqcData(x)
              plotNContent(x, usePlotly, labels, pwfCols, warn, fail, ...)
          }
)
#' @aliases plotNContent,FastqcFileList
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, ...){
              x <- getFastqcData(x)
              plotNContent(x, usePlotly, labels, pwfCols, warn, fail, ...)
          }
)
#' @aliases plotNContent,FastqcData
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, 
                   ..., lineCol = "red"){

              # Get the NContent
              df <- Per_base_N_content(x)

              # Handle empty/missing modules
              if (!length(df)) {
                  #stop("N Content Module")
                  nPlot <- emptyPlot("No N Content Module Detected")

                  if(usePlotly) nPlot <- ggplotly(nPlot, tooltip = "")
                  return(nPlot)
              }
              if (sum(df[["N-Count"]]) == 0) {
                  #stop("No N content was detected by FastQC")
                  nPlot <- ngsReports:::emptyPlot("No N Content in Sequences")

                  if(usePlotly) nPlot <- ggplotly(nPlot, tooltip = "")
                  return(nPlot)
              }

              colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

              # Sort out the colours
              if (missing(pwfCols)) pwfCols <- ngsReports::pwf
              stopifnot(isValidPwf(pwfCols))
              pwfCols <- setAlpha(pwfCols, 0.2)

              labels <- makeLabels(df, labels, ...)
              df$Filename <- labels[df$Filename]
              df$Base <- factor(df$Base, levels = unique(df$Base))
              df$xValue <- as.integer(df$Base)

              # Setup the BG colours
              rects <- tibble::tibble(xmin = 0,
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

              yLab <- "N Content (%)"
              nPlot <- ggplot(df) +
                  geom_rect(data = rects,
                            aes_string(xmin = "xmin", xmax = "xmax",
                                       ymin = "ymin", ymax = "ymax",
                                       fill = "Status")) +
                  geom_line(aes_string(x = "xValue", y = "Percentage"),
                            colour = lineCol) +
                  geom_point(aes_string(x = "xValue", y = "Percentage",
                                        group = "Base"),
                             size = 0, colour = rgb(0, 0, 0, 0)) +
                  scale_fill_manual(values = getColours(pwfCols)) +
                  scale_x_continuous(breaks = unique(df$xValue),
                                     labels = levels(df$Base), expand = c(0,0)) +
                  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                  facet_wrap(~Filename) +
                  labs(x = "Position in Read",
                       y = yLab) +
                  guides(fill = FALSE) +
                  theme_bw() +
                  theme(axis.text.x = element_text(
                      angle = 90, hjust = 1, vjust = 0.5
                      ))

              # Add the basic customisations
              if (!is.null(userTheme)) nPlot <- nPlot + userTheme

              if (usePlotly){
                  nPlot <- nPlot +
                      xlab("") +
                      theme(legend.position = "none")
                  nPlot <- suppressMessages(plotly::ggplotly(nPlot))
                      # ))

                  nPlot <- suppressMessages(
                      suppressWarnings(
                          plotly::subplot(plotly::plotly_empty(),
                                          nPlot, widths = c(0.14,0.86))
                      ))
                  nPlot <- plotly::layout(nPlot, yaxis2 = list(title = yLab))


                  # Set the hoverinfo for bg rectangles to the vertices only,
                  # This will effectively hide them
                  nPlot$x$data[[1]]$hoverinfo <- "none"
                  nPlot$x$data[[2]]$hoverinfo <- "none"
                  nPlot$x$data[[3]]$hoverinfo <- "none"
                  nPlot$x$data[[4]]$hoverinfo <- "none"
                  nPlot$x$data[[5]]$hoverinfo <- "none"
                  # Hide the xValue parameter to make it look nicer
                  nPlot$x$data[[6]]$text <- gsub("(.+)(xValue.+)(Percentage.+)",
                                                 "\\1\\3", nPlot$x$data[[6]]$text)

              }

              nPlot

          }
)
#' @aliases plotNContent,FastqcDataList
#' @rdname plotNContent-methods
#' @export
setMethod("plotNContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 20, 
                   cluster = FALSE, dendrogram = FALSE, ...){
              # Get the NContent
              df <- Per_base_N_content(x)


              if (!length(df)) {
                  nPlot <- emptyPlot("No N Content Module Detected")
                  if(usePlotly) acPlot <- ggplotly(nPlot, tooltip = "")
                  return(nPlot)
              }

              if (sum(df$`N-Count`) == 0) {
                  nPlot <- ngsReports:::emptyPlot("No N Content in Sequences")
                  if(usePlotly) acPlot <- ggplotly(nPlot, tooltip = "")
                  return(nPlot)
              }

              colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

              # Sort out the colours
              if (missing(pwfCols)) pwfCols <- ngsReports::pwf
              stopifnot(isValidPwf(pwfCols))

              # Get the labels organised
              labels <- makeLabels(df, labels, ...)

              ## fill bins up to the max sequence length
              df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))
              df <- lapply(split(df, f = df$Filename),
                           function(x){
                               Longest_sequence <- max(x$Start)
                               dfFill <- data.frame(Start = seq_len(Longest_sequence))
                               x <- dplyr::right_join(x, dfFill, by = "Start")
                               x <- na.locf(x)
                           })
              df <- dplyr::bind_rows(df)

              if (dendrogram && !cluster){
                  message("cluster will be set to TRUE when dendrogram = TRUE")
                  cluster <- TRUE
              }

              # Now define the order for a dendrogram if required
              key <- names(labels)
              if (cluster){
                  clusterDend <- setClusters(df = df[c("Filename", "Start", "Percentage")],
                                             rowVal = "Filename",
                                             colVal = "Start",
                                             value = "Percentage")
                  key <- labels(clusterDend)
              }
              # Now set everything as factors
              df$Filename <- factor(labels[df$Filename],
                                    levels = labels[key])

              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

              xLab <- "Position in Read (bp)"

              nPlot <- ggplot(df, aes_string("Start", "Filename",
                                             fill = "Percentage", label = "Base")) +
                  geom_tile() +
                  scale_fill_pwf(df$Percentage, pwfCols,
                                 breaks = c(0, warn, fail, 101),
                                 passLow = TRUE, na.value = "white") +
                  scale_x_continuous(expand = c(0,0),
                                     breaks = unique(df$x),
                                     labels = unique(df$Base)) +
                  scale_y_discrete(expand = c(0, 0)) +
                  labs(x = xLab,
                       y = "Filename",
                       fill = "%N") +
                  theme_bw() +
                  theme(panel.background = element_blank(),
                        axis.text.x = element_text(angle = 90,
                                                   hjust = 1,
                                                   vjust = 0.5))


              # Add the basic customisations
              if (!is.null(userTheme)) nPlot <- nPlot + userTheme

              if (usePlotly){
                  # Reset the status using current values
                  status <- dplyr::summarise_at(dplyr::group_by(df, Filename),
                                                vars("Percentage"),
                                                funs(Percentage = max), na.rm = TRUE)
                  status$Status <- cut(status$Percentage,
                                       breaks = c(0, warn, fail, 101),
                                       include.lowest = TRUE,
                                       labels = c("PASS", "WARN", "FAIL"))

                  # Form the sideBar for each adapter
                  sideBar <- makeSidebar(status, key, pwfCols = pwfCols)

                  if (!is.null(userTheme)) nPlot <- nPlot + userTheme

                  if (dendrogram){
                      dx <- ggdendro::dendro_data(clusterDend)
                      dendro <- ggdend(dx$segments) +
                          coord_flip() +
                          scale_y_reverse(expand = c(0, 0)) +
                          scale_x_continuous(expand = c(0, 0.5))
                      dendro <- plotly::ggplotly(dendro, tooltip = NULL)
                      nPlot <- nPlot + ylab("")
                  }
                  else{
                      dendro <- plotly::plotly_empty()
                  }

                  # Customise for plotly
                  nPlot <- nPlot +
                      theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank())
                  nPlot <- plotly::ggplotly(nPlot,
                                            tooltip = c("y", "fill", "label"))
                  # Now make the three frame plot with option dendrogram
                  # and the sideBar + main plot
                  nPlot <- suppressWarnings(
                      suppressMessages(
                          plotly::subplot(dendro, sideBar, nPlot,
                                          widths = c(0.1,0.08,0.82),
                                          margin = 0.001,
                                          shareY = TRUE)
                      ))
                  nPlot <- plotly::layout(nPlot,
                                          xaxis3 = list(title = xLab),
                                          margin = list(b = 45))

              }

              nPlot

          }
)
