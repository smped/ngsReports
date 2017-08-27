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
#' @param adapterType A regular expression matching the adapter(s) to be plotted
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param warn,fail The default values for warn and fail are 5 and 10 respectively (i.e. precentages)
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
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
#' plotAdapterContent(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fileName(fdl))
#' plotAdapterContent(fdl[r1])
#'
#' # Plot just the Universal Adapter
#' # and change the y-axis using ggplot2::scale_y_continuous
#' library(ggplot2)
#' plotAdapterContent(fdl, adapterType ="Universal", plotType = "line") +
#' scale_y_continuous()
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
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
#' @importFrom dplyr funs
#'
#' @name plotAdapterContent
#' @rdname plotAdapterContent-methods
#' @export
setGeneric("plotAdapterContent",function(x, usePlotly = FALSE, ...){standardGeneric("plotAdapterContent")})
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
          function(x, usePlotly = FALSE, pwfCols, warn = 5, fail = 10, ...){

            df <- Adapter_Content(x)

            # Sort out the colours & pass/warn/fail breaks
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))
            stopifnot(is.numeric(c(warn, fail)))
            stopifnot(all(fail < 100, warn < fail,  warn > 0))

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Change to long form and remove the _ symbols between words
            df <- reshape2::melt(df, id.vars = c("Filename", "Position"),
                                 value.name = "Percent", variable.name = "Type")
            df <- dplyr::mutate(df, Type = gsub("_", " ", Type))

            # Set the positions as a factor
            df$Position <- gsub("([0-9]*)-.+", "\\1", df$Position)
            df$Position <- as.numeric(df$Position)
            # Round the percent for nicer plotting
            df$Percent <- round(df$Percent, 2)

            # Add transparency to background colours & define the rectangles
            pwfCols <- setAlpha(pwfCols, 0.2)
            x <- list(min = min(df$Position), max = max(df$Position))
            rects <- data_frame(xmin = 0,
                                xmax = max(df$Position),
                                ymin = c(0, warn, fail),
                                ymax = c(warn, fail, 100),
                                Status = c("PASS", "WARN", "FAIL"))

            # Create the basic plot
            acPlot <- ggplot(df) +
              geom_rect(data = rects,
                        aes_string(xmin = "xmin", xmax = "xmax",
                                   ymin = "ymin", ymax = "ymax", fill = "Status")) +
              geom_line(aes_string(x = "Position", y = "Percent",colour = "Type")) +
              scale_y_continuous(limits = c(0, 100), expand =c(0, 0)) +
              scale_x_continuous(expand = c(0, 0)) +
              scale_fill_manual(values = getColours(pwfCols)) +
              scale_colour_discrete() +
              facet_wrap(~Filename, ncol = 1) +
              labs(x = "Position",
                   y = "Percent (%)") +
              guides(fill = FALSE) +
              theme_bw() +
              theme(legend.position = c(1, 1),
                    legend.justification = c(1, 1),
                    legend.title = element_blank())

            # Add the basic customisations
            if (!is.null(userTheme)) acPlot <- acPlot + userTheme
            # Make interacive if required
            if (usePlotly){
              acPlot <- suppressMessages(
                plotly::ggplotly(acPlot + theme(legend.position = "none"),
                                 hoverinfo = c("x", "y", "colour"))
              )
              # Set the hoverinfo for bg rectangles to the vertices only,
              # This will effectively hide them
              acPlot$x$data[[1]]$hoveron <- "points"
              acPlot$x$data[[2]]$hoveron <- "points"
              acPlot$x$data[[3]]$hoveron <- "points"
            }

            acPlot

          }
)
#' @aliases plotAdapterContent,FastqcDataList
#' @rdname plotAdapterContent-methods
#' @export
setMethod("plotAdapterContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, plotType = "heatmap", labels, adapterType,
                   pwfCols, warn = 5, fail = 10, ...){

            df <- Adapter_Content(x)

            # Sort out the colours & pass/warn/fail breaks
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))
            stopifnot(is.numeric(c(warn, fail)))
            stopifnot(all(fail < 100, warn < fail,  warn > 0))

            # Check for valid plotType
            stopifnot(plotType %in% c("line", "heatmap"))

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)),
                                  names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Set the positions as a factor
            df$Position <- factor(df$Position, levels = unique(df$Position))

            # Change to long form and remove the _ symbols between words
            df <- reshape2::melt(df, id.vars = c("Filename", "Position"),
                                 value.name = "Percent", variable.name = "Type")
            df <- dplyr::mutate(df, Type = gsub("_", " ", Type))

            # Restrict to a given adapterType if requested
            if (!missing(adapterType)) {
              keep <- grep(pattern = adapterType, x = df$Type)
              df <- df[keep,]
              if(nrow(df) == 0) stop("No adapters matching the supplied type were found")
            }

            # Decide if plots should be dropped to avoid issues with plotly
            # and heatmaps which contain all zeroes
            maxAdapters <- dplyr::group_by(df, Type)
            maxAdapters <- dplyr::summarise_at(maxAdapters, vars("Percent"), funs(Percent = max))
            keepType <- maxAdapters$Type[maxAdapters$Percent > 0]
            if (length(keepType) != nrow(maxAdapters)){
              dropType <- maxAdapters$Type[!maxAdapters$Percent > 0]
              message("Adapters with all zero values found. Removing plots for:\n",
                      paste(dropType, collapse = "\n"))
            }
            if (length(keepType) == 0) {
              message("Adapter content is all zero. No plots are possible")
              return(NULL)
            }
            df <- dplyr::filter(df, Type %in% keepType)
            df <- droplevels(df)

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            breaks <- c(0, warn, fail, 100)

            if (plotType == "heatmap"){

              # Reverse the factor levels for a better looking plot
              df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

              if (!usePlotly){
                acPlot <- ggplot(df, aes_string(x = "Position", y = "Filename", fill = "Percent")) +
                  scale_y_discrete(labels = labels) +
                  geom_tile() +
                  facet_wrap(~Type, ncol = 1) +
                  scale_fill_pwf(df$Percent, pwf, breaks = breaks)
              }
              else{
                # Now generate plots for each adapter
                splitDf <- split(df, f = df$Type)
                acPlotList <- lapply(splitDf, function(x){

                  # Reset the status using adapter specific values
                  status <- dplyr::summarise_at(dplyr::group_by(x, Filename), vars("Percent"), funs(Percent = max))
                  status$Status <- cut(status$Percent, breaks = breaks, include.lowest = TRUE,
                                       labels = c("PASS", "WARN", "FAIL"))
                  status$x <- 1

                  # Form the sideBar for each adapter
                  sideBar <- ggplot(status, aes_string("x", "Filename")) +
                    geom_tile(aes_string(fill = "Status")) +
                    scale_fill_manual(values = getColours(pwfCols)) +
                    scale_y_discrete(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0, 0)) +
                    theme(panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          legend.position="none",
                          axis.title=element_blank(),
                          axis.text=element_blank(),
                          axis.ticks=element_blank())
                  sideBar <- suppressMessages(ggplotly(sideBar, tooltip = c("Status", "Filename")))

                  # Make the individual heatmap
                  p <- ggplot(x, aes_string(x = "Position", y = "Filename", fill = "Percent")) +
                    scale_y_discrete(labels = labels) +
                    geom_tile() +
                    scale_fill_pwf(x$Percent, pwfCols, breaks =c(0, warn, fail, 100)) +
                    facet_wrap(~Type, ncol = 1) +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())

                  # Return each plotly object into a list
                  suppressMessages(
                    plotly::subplot(sideBar, p, widths = c(0.1,0.9), margin = 0.01, shareY = TRUE)
                  )
                })
                acPlot <- suppressMessages(
                  plotly::subplot(acPlotList, margin = c(0,0,0.06,0.06), titleX = TRUE,
                                  nrows = length(splitDf))
                )
              }
            }

            if (plotType == "line") {

              # Define position as an integer
              df$Position <- gsub("([0-9]*)-.+", "\\1", df$Position)
              df$Position <- as.numeric(df$Position)

              # Set the transparency & position of bg rectangles
              pwfCols <- setAlpha(pwfCols, 0.2)
              rects <- data_frame(xmin = 0,
                                  xmax = max(df$Position),
                                  ymin = c(0, warn, fail),
                                  ymax = c(warn, fail, 100),
                                  Status = c("PASS", "WARN", "FAIL"))

              # Create the basic plot
              acPlot <- ggplot(df) +
                geom_rect(data = rects,
                          aes_string(xmin = "xmin", xmax = "xmax",
                                     ymin = "ymin", ymax = "ymax", fill = "Status")) +
                geom_line(aes_string(x = "Position", y = "Percent", colour = "Filename")) +
                scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_colour_discrete(labels = labels) +
                scale_fill_manual(values = getColours(pwfCols)) +
                guides(fill = FALSE) +
                labs(x = "Position in read (bp)",
                     y = "Percent (%)") +
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
                acPlot$x$data[[1]]$hoveron <- "points"
                acPlot$x$data[[2]]$hoveron <- "points"
                acPlot$x$data[[3]]$hoveron <- "points"
              }

            }
            acPlot
          }
)
