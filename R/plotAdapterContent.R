#' @title Draw an Adapter Content Plot
#'
#' @description Draw an Adapter Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the Adapter_Content from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' When \code{x} is a single FastqcFile, or FastqcData object line plots will always
#' be drawn for all adapters.
#' Otherwise, users can select line plots or heatmaps.
#' When plotting more than one fastqc file, any undetected adapters will not be shown.
#'
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical}. Output as ggplot2 (default) or plotly object.
#' @param adapterType A regular expression matching the adapter(s) to be plotted.
#' To plot all adapters summed, specify \code{adapterType = "Total"}.
#' This is the default behaviour.
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param warn,fail The default values for warn and fail are 5 and 10 respectively (i.e. precentages)
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
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
#' plotAdapterContent(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fileName(fdl))
#' plotAdapterContent(fdl[r1])
#'
#' # Plot just the Universal Adapter
#' # and change the y-axis using ggplot2::scale_y_continuous
#' plotAdapterContent(fdl, adapterType ="Universal", plotType = "line") +
#' facet_wrap(~Filename) +
#' guides(colour = FALSE)
#'
#' @import ggplot2
#' @import tibble
#' @importFrom plotly plotly_empty ggplotly
#' @importFrom stats hclust dist
#' @importFrom zoo na.locf
#'
#' @name plotAdapterContent
#' @rdname plotAdapterContent-methods
#' @export
setGeneric("plotAdapterContent",function(x, usePlotly = FALSE, ...){
    standardGeneric("plotAdapterContent")
})
#' @aliases plotAdapterContent,character
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "character",
          function(x, usePlotly = FALSE, ...){
              x <- getFastqcData(x)
              plotAdapterContent(x, usePlotly,...)
          }
)
#' @aliases plotAdapterContent,FastqcFile
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
              x <- getFastqcData(x)
              plotAdapterContent(x, usePlotly,...)
          }
)
#' @aliases plotAdapterContent,FastqcFileList
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
              x <- getFastqcData(x)
              plotAdapterContent(x, usePlotly,...)
          }
)
#' @aliases plotAdapterContent,FastqcData
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, pwfCols, warn = 5, fail = 10, ...){

              df <- Adapter_Content(x)

              if (!length(df)) {

                  acPlot <- emptyPlot("No Adapter Content Module Detected")
                  if(usePlotly) acPlot <- ggplotly(acPlot, tooltip = "")

                  return(acPlot)
              }

              valueCols <- setdiff(colnames(df), c("Filename", "Position"))
              if (sum(df[valueCols]) == 0) {

                  acPlot <- emptyPlot("No Adapter Content Found in Sequences")
                  if(usePlotly) acPlot <- ggplotly(acPlot, tooltip = "")

                  return(acPlot)
              }

              # Set any labels
              labels <- setLabels(df, labels, ...)
              df$Filename <- labels[df$Filename]

              # Sort out the colours & pass/warn/fail breaks
              if (missing(pwfCols)) pwfCols <- ngsReports::pwf
              stopifnot(isValidPwf(pwfCols))
              stopifnot(is.numeric(c(warn, fail)))
              stopifnot(all(fail < 100, warn < fail,  warn > 0))

              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme,
                                                             dotArgs[keepArgs])

              # Change to long form and remove the _ symbols between words
              df <- tidyr::gather(df, key = "Type", value = "Percent",
                                  tidyselect::one_of(valueCols))
              df$Type <- gsub("_", " ", df$Type)

              # Set the positions as a factor
              df$Position <- gsub("([0-9]*)-.+", "\\1", df$Position)
              df$Position <- as.numeric(df$Position)
              # Round the percent for nicer plotting with plotly
              df$Percent <- round(df$Percent, 2)

              # Add transparency to background colours & define the rectangles
              pwfCols <- setAlpha(pwfCols, 0.2)
              rng <- structure(range(df$Position), names = c("min", "max"))
              rects <- tibble(xmin = 0,
                                      xmax = rng["max"],
                                      ymin = c(0, warn, fail),
                                      ymax = c(warn, fail, 100),
                                      Status = c("PASS", "WARN", "FAIL"))

              # Create the basic plot
              xlab <- "Position in read (bp)"
              ylab <- "Percent (%)"
              acPlot <- ggplot(df) +
                  geom_rect(data = rects,
                            aes_string(xmin = "xmin", xmax = "xmax",
                                       ymin = "ymin", ymax = "ymax",
                                       fill = "Status")) +
                  geom_line(aes_string(x = "Position", y = "Percent",
                                       colour = "Type")) +
                  scale_y_continuous(limits = c(0, 100), expand =c(0, 0)) +
                  scale_x_continuous(expand = c(0, 0)) +
                  scale_fill_manual(values = getColours(pwfCols)) +
                  scale_colour_discrete() +
                  facet_wrap(~Filename, ncol = 1) +
                  labs(x = xlab,
                       y = ylab) +
                  guides(fill = FALSE) +
                  theme_bw() +
                  theme(legend.position = c(1, 1),
                        legend.justification = c(1, 1),
                        legend.title = element_blank())

              # Add the basic customisations
              if (!is.null(userTheme)) acPlot <- acPlot + userTheme
              # Make interacive if required
              if (usePlotly){
                  acPlot <- acPlot + theme(legend.position = "none")
                  acPlot <- suppressWarnings(
                      suppressMessages(
                          plotly::subplot(
                              plotly::plotly_empty(),
                              acPlot,
                              widths = c(0.14,0.86)))
                  )
                  acPlot <- plotly::layout(acPlot,
                                           xaxis2 = list(title = xlab),
                                           yaxis2 = list(title = ylab))
                  # Set the hoverinfo for bg rectangles to the vertices only,
                  # This will effectively hide them
                  acPlot$x$data[[1]]$hoverinfo <- "none"
                  acPlot$x$data[[2]]$hoverinfo <- "none"
                  acPlot$x$data[[3]]$hoverinfo <- "none"
              }

              acPlot

          }
)
#' @aliases plotAdapterContent,FastqcDataList
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, plotType = c("heatmap", "line"), labels,
                   adapterType = "Total",
                   pwfCols, warn = 5, fail = 10,
                   cluster = FALSE, dendrogram = FALSE, ...){

              df <- Adapter_Content(x)

              # Check the two conditions for an empty plot
              if (!length(df)) {
                  acPlot <- emptyPlot("No Adapter Content Module Detected")
                  if(usePlotly) acPlot <- ggplotly(acPlot, tooltip = "")
                  return(acPlot)
              }

              valueCols <- setdiff(colnames(df), c("Filename", "Position"))
              if (sum(df[valueCols]) == 0) {
                  acPlot <- emptyPlot("No Adapter Content in Sequences")
                  if(usePlotly) acPlot <- ggplotly(acPlot, tooltip = "")
                  return(acPlot)
              }

              # Sort out the colours & pass/warn/fail breaks
              if (missing(pwfCols)) pwfCols <- ngsReports::pwf
              stopifnot(isValidPwf(pwfCols))
              stopifnot(is.numeric(c(warn, fail)))
              stopifnot(all(fail < 100, warn < fail,  warn > 0))
              breaks <- c(0, warn, fail, 100)

              # Check for valid plotType & labels
              plotType <- match.arg(plotType[1], c("line", "heatmap"))
              labels <- setLabels(df, labels, ...)

              # Change to long form
              df <- tidyr::gather(df, key = "Type", value = "Percent",
                                  tidyselect::one_of(valueCols))
              # Set the position as a factor
              df$Position <- factor(df$Position, levels = unique(df$Position))

              # Get the adapter type and summarise all if total is requested
              adapterType <- grep(adapterType[1], c("Total_Adapter_Content", valueCols),
                                  value = TRUE)
              if (length(adapterType) != 1) stop("Could not determine adapter type")
              if (adapterType == "Total_Adapter_Content") {
                  # Sum the adapters by filename& position
                  df <- dplyr::group_by(df, Filename, Position)
                  df <- dplyr::summarise_at(df,
                                            dplyr::vars("Percent"),
                                            dplyr::funs(Percent = sum),
                                            na.rm = TRUE)
                  df <- dplyr::ungroup(df)
              }
              else{
                  df <- dplyr::filter(df, Type == adapterType)
              }

              # Remove the underscores from the adapterType for prettier output
              adapterType <- gsub("_", " ", adapterType)

              # If no adapter content is found, output a simple message
              if (max(df$Percent) == 0){
                  msg <- paste("No", adapterType, "found")
                  return(emptyPlot(msg))
              }

              # Now just keep the three required columns
              df <- df[c("Filename", "Position", "Percent")]
              df <- droplevels(df)
              df$Percent <- round(df$Percent, 4)

              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

              # Set the axis label for either plotType
              xlab <- "Position in read (bp)"

              if (plotType == "heatmap"){

                  ylab <- ifelse(dendrogram, "", "Filename")

                  # Get the longest sequence
                  df$Start <- as.integer(
                      gsub("([0-9]*)-[0-9]*", "\\1", as.character(df$Position))
                  )
                  df <- df[c("Filename", "Start", "Percent", "Position")]

                  # Adjust the data for files with varying read lengths
                  # This will fill NA values with the previous values
                  df <- lapply(split(df, f = df$Filename), function(x){
                      Longest_sequence <- gsub(".*-([0-9]*)", "\\1",
                                               as.character(x$Position))
                      Longest_sequence <- max(as.integer(Longest_sequence))
                      dfFill <- data.frame(Start = seq_len(Longest_sequence))
                      x <- dplyr::right_join(x, dfFill, by = "Start")
                      na.locf(x)
                  })
                  df <- dplyr::bind_rows(df)

                  # Use the start position for plotting instead of Position
                  plotCols <- c("Filename", "Start", "Percent")
                  df <- df[plotCols]

                  # Arrange by row if clustering
                  # Just use the default order as the key if not clustering
                  # Always turn clustering on if dendrogram = TRUE
                  if (dendrogram && !cluster){
                      message("cluster will be set to TRUE when dendrogram = TRUE")
                      cluster <- TRUE
                  }
                  key <- names(labels)
                  if(cluster){
                      clusterDend <- setClusters(df, "Filename", "Start", "Percent")
                      key <- labels(clusterDend)
                  }

                  # Set the factor levels for the y-axis
                  df$Filename <- factor(labels[df$Filename], levels = labels[key])
                  df$Type <- adapterType

                  # Make the heatmap
                  acPlot <- ggplot(df,
                                   aes_string(x = "Start", y = "Filename",
                                              fill = "Percent", type = "Type")) +
                      geom_tile() +
                      ggtitle(adapterType) +
                      labs(x = xlab) +
                      scale_x_continuous(expand = c(0,0)) +
                      scale_y_discrete(expand = c(0, 0)) +
                      scale_fill_pwf(df$Percent, pwf, breaks = breaks,
                                     na.value = "white") +
                      theme_bw() +
                      theme(plot.title = element_text(hjust = 0.5))

                  if (usePlotly){

                      # Reset the PWF status using current values
                      # This needs to be recalculated when using Total AC
                      status <- dplyr::summarise_at(
                          dplyr::group_by(df, Filename),
                          dplyr::vars("Percent"),
                          dplyr::funs(Percent = max),
                          na.rm = TRUE)
                      status$Status <- cut(status$Percent, breaks = breaks,
                                           include.lowest = TRUE,
                                           labels = c("PASS", "WARN", "FAIL"))

                      # Form the sideBar for each adapter
                      sideBar <- makeSidebar(status, key, pwfCols = pwfCols)

                      # Customise for plotly
                      acPlot <- acPlot +
                          ggtitle(NULL) +
                          theme(panel.background = element_blank(),
                                axis.text.x = element_text(angle = 90,
                                                           hjust = 1,
                                                           vjust = 0.5),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                legend.position = "none")
                      if (!is.null(userTheme)) acPlot <- acPlot + userTheme

                      # plot dendro
                      if (dendrogram){
                          dx <- ggdendro::dendro_data(clusterDend)
                          dendro <- ggdend(dx$segments) +
                              coord_flip() +
                              scale_y_reverse(expand = c(0, 0)) +
                              scale_x_continuous(expand = c(0, 0.5))
                          dendro <- plotly::ggplotly(dendro, tooltip = NULL)
                      }
                      else {
                          dendro <- plotly::plotly_empty()
                      }

                      acPlot <- suppressWarnings(
                          suppressMessages(
                              plotly::subplot(dendro,
                                              sideBar,
                                              acPlot,
                                              widths = c(0.1,0.08,0.82),
                                              margin = 0.001,
                                              shareY = TRUE)
                          ))
                      acPlot <- plotly::layout(acPlot,
                                               annotations = list(
                                                   text = ylab,
                                                   textangle = -90,
                                                   showarrow = FALSE
                                               ),
                                               xaxis3 = list(title = xlab),
                                               margin = list(b = 52))

                  }
              }

              if (plotType == "line") {

                  ylab <- "Percent (%)"

                  # Define position as an integer just taking the first
                  # value in the Position field
                  df$Start <- gsub("([0-9]*)-.+", "\\1",
                                   as.character(df$Position))
                  df$Start <- as.integer(df$Start)
                  df$Position <- as.numeric(df$Position)
                  df$Type <- adapterType

                  # Set the transparency & position of bg rectangles
                  pwfCols <- setAlpha(pwfCols, 0.2)
                  rects <- tibble(xmin = 0,
                                          xmax = max(df$Start),
                                          ymin = c(0, warn, fail),
                                          ymax = c(warn, fail, 100),
                                          Status = c("PASS", "WARN", "FAIL"))

                  # Create the basic plot
                  acPlot <- ggplot(df) +
                      geom_rect(data = rects,
                                aes_string(xmin = "xmin", xmax = "xmax",
                                           ymin = "ymin", ymax = "ymax",
                                           fill = "Status")) +
                      geom_line(aes_string(x = "Start", y = "Percent",
                                           colour = "Filename")) +
                      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_colour_discrete(labels = labels) +
                      scale_fill_manual(values = getColours(pwfCols)) +
                      guides(fill = FALSE) +
                      labs(x = xlab, y = ylab) +
                      facet_wrap(~Type, ncol = 1) +
                      theme_bw()
                  if (!is.null(userTheme)) acPlot <- acPlot + userTheme

                  # And draw the plot
                  if (usePlotly){
                      acPlot <- acPlot + theme(legend.position = "none")
                      acPlot <- suppressMessages(
                          plotly::ggplotly(acPlot,
                                           hoverinfo = c("x", "y", "colour"))
                      )

                      # Set the hoverinfo for bg rectangles to the vertices only,
                      # This will effectively hide them
                      acPlot$x$data <- lapply(acPlot$x$data, function(x){
                          # If there is a name component & it contains PASS/WARN/FAIL
                          # set the hoverinfo to none
                          if ("name" %in% names(x)){
                              if (grepl("(PASS|WARN|FAIL)", x$name)){
                                  x$hoverinfo <- "none"
                              }
                          }
                          x
                      })
                  }
              }
              acPlot
          }
)
