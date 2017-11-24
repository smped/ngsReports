#' @title Plot the Base Qualities for each file
#'
#' @description Plot the Base Qualities for each file as separate plots
#'
#' @details This replicates the \code{Per base sequence quality} plots from FASTQC,
#' using facets to plce them all in a single ggplot2 object.
#'
#' For large datasets, subsetting by R1 or R2 reads may be helpful
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param nc \code{numeric}. The number of columns to create in the plot layout
#' @param warn,fail The default values for warn and fail are 30 and 20 respectively (i.e. percentages)
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param plotType \code{character} Can be either \code{"boxplot"} or \code{"heatmap"}
#' @param plotValue \code{character} Type of quality data to be presented "Mean" or "Median"
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A standard ggplot2 object
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
#' # The default plot for multiple libraries is a heatmap
#' plotBaseQualities(fdl)
#'
#' # The default plot for a single library is the standard boxplot
#' plotBaseQualities(fdl[[1]])
#'
#' @import ggplot2
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom stats na.omit
#' @importFrom stats hclust dist
#'
#' @name plotBaseQualities
#' @rdname plotBaseQualities-methods
#' @export
setGeneric("plotBaseQualities",function(x, usePlotly = FALSE, ...){standardGeneric("plotBaseQualities")})
#' @aliases plotBaseQualities,character
#' @rdname plotBaseQualities-methods
#' @export
setMethod("plotBaseQualities", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotBaseQualities(x, usePlotly,...)
          }
)
#' @aliases plotBaseQualities,FastqcFile
#' @rdname plotBaseQualities-methods
#' @export
setMethod("plotBaseQualities", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotBaseQualities(x, usePlotly,...)
          }
)
#' @aliases plotBaseQualities,FastqcFileList
#' @rdname plotBaseQualities-methods
#' @export
setMethod("plotBaseQualities", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotBaseQualities(x, usePlotly,...)
          }
)
#' @aliases plotBaseQualities,FastqcData
#' @rdname plotBaseQualities-methods
#' @export
setMethod("plotBaseQualities", signature = "FastqcData",
          function(x, usePlotly = FALSE, pwfCols, warn = 25, fail = 20, ...){

            # Get the data
            df <- tryCatch(Per_base_sequence_quality(x))
            df$Base <- factor(df$Base, levels = unique(df$Base))
            df$x <- as.integer(df$Base)
            df$xmin <- df$x - 0.4
            df$xmax <- df$x + 0.4

            # Sort out the colours
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))
            pwfCols <- setAlpha(pwfCols, 0.2)

            # Get any theme arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            # Set the limits & rectangles
            ylim <- c(0, max(df$`90th_Percentile`) + 1)
            expand_x <- round(0.015*(max(df$x) - min(df$x)), 1)
            rects <- dplyr::data_frame(xmin = min(df$x) - expand_x,
                                       xmax = max(df$x) + expand_x,
                                       ymin = c(0, fail, warn),
                                       ymax = c(fail, warn, max(ylim)),
                                       Status = c("FAIL", "WARN", "PASS"))

            # Get the Illumina encoding
            enc <- Basic_Statistics(x)$Encoding[1]
            enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)

            # Generate the basic plot
            qualPlot <- ggplot(df) +
              geom_rect(data = rects,
                        aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax",
                                   fill = "Status")) +
              geom_rect(aes_string(xmin = "xmin", xmax = "xmax",
                                   ymin = "Lower_Quartile", ymax = "Upper_Quartile"),
                        fill ="yellow", colour = "black") +
              geom_segment(aes_string(x = "xmin", xend = "xmax", y = "Median", yend = "Median"),
                           colour = "red") +
              geom_linerange(aes_string(x = "x", ymin = "`10th_Percentile`", ymax = "Lower_Quartile")) +
              geom_linerange(aes_string(x = "x", ymin = "Upper_Quartile", ymax = "`90th_Percentile`")) +
              geom_line(aes_string(x = "x", y = "Mean"), colour = "blue") +
              scale_fill_manual(values = getColours(pwfCols)) +
              scale_x_continuous(breaks = df$x, labels = levels(df$Base), expand = c(0, 0)) +
              scale_y_continuous(limits = ylim, expand = c(0, 0)) +
              xlab("Position in read (bp)") +
              ylab(paste0("Quality Scores (", enc, " encoding)")) +
              facet_wrap(~Filename, ncol = 1) +
              guides(fill = FALSE) +
              theme_bw() +
              theme(panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
            if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

            if(usePlotly){
              qualPlot <- qualPlot +
                xlab("") +
                theme(legend.position = "none")
              qualPlot <- suppressMessages(
                plotly::ggplotly(qualPlot, hoverinfo = c("Base", "Mean", "Median",
                                                         "Upper_Quartile", "Lower_Quartile",
                                                         "`10th_Percentile`","`90th_Percentile`"))
              )
            }

            qualPlot
          }
)
#' @aliases plotBaseQualities,FastqcDataList
#' @rdname plotBaseQualities-methods
#' @export
setMethod("plotBaseQualities", signature = "FastqcDataList",
          function(x,  usePlotly = FALSE,  plotType = c("heatmap", "boxplot"), plotValue = "Mean",
                   cluster = FALSE, labels, dendrogram = FALSE,
                   nc = 2, pwfCols, warn = 25, fail = 20, ...){

            # Get the data
            df <- Per_base_sequence_quality(x)

            # Sort out the colours
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf
            stopifnot(isValidPwf(pwfCols))

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Get the Illumina encoding
            enc <- Basic_Statistics(x)$Encoding[1]
            enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)

            plotType <- match.arg(plotType)

            if (plotType == "boxplot"){

              # Sort out the colours
              pwfCols <- setAlpha(pwfCols, 0.2)

              # Setup the boxes & the x-axis
              df$Base <- factor(df$Base, levels = unique(df$Base))
              df$x <- as.integer(df$Base)
              df$xmin <- df$x - 0.4
              df$xmax <- df$x + 0.4

              # Set the limits & rectangles
              ylim <- c(0, max(df$`90th_Percentile`) + 1)
              expand_x <- round(0.015*(max(df$x) - min(df$x)), 1)
              rects <- dplyr::data_frame(xmin = min(df$x) - expand_x,
                                         xmax = max(df$x) + expand_x,
                                         ymin = c(0, fail, warn),
                                         ymax = c(fail, warn, max(ylim)),
                                         Status = c("FAIL", "WARN", "PASS"))

              # Get any theme arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

              # Generate the basic plot
              df$Filename <- labels[df$Filename]
              qualPlot <- ggplot(df) +
                geom_rect(data = rects,
                          aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax",
                                     fill = "Status")) +
                geom_rect(aes_string(xmin = "xmin", xmax = "xmax",
                                     ymin = "Lower_Quartile", ymax = "Upper_Quartile"),
                          fill ="yellow", colour = "black") +
                geom_segment(aes_string(x = "xmin", xend = "xmax", y = "Median", yend = "Median"),
                             colour = "red") +
                geom_linerange(aes_string(x = "x", ymin = "`10th_Percentile`", ymax = "Lower_Quartile")) +
                geom_linerange(aes_string(x = "x", ymin = "Upper_Quartile", ymax = "`90th_Percentile`")) +
                geom_line(aes_string(x = "x", y = "Mean"), colour = "blue") +
                scale_fill_manual(values = getColours(pwfCols)) +
                scale_x_continuous(breaks = unique(df$x), labels = levels(df$Base), expand = c(0, 0)) +
                scale_y_continuous(limits = ylim, expand = c(0, 0)) +
                xlab("Position in read (bp)") +
                ylab(paste0("Quality Scores (", enc, " encoding)")) +
                facet_wrap(~Filename, ncol = nc) +
                guides(fill = FALSE) +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
              if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

              # Make interactive if required
              if(usePlotly){
                qualPlot <- qualPlot +
                  xlab("") +
                  theme(legend.position = "none")
                qualPlot <- suppressMessages(
                  plotly::ggplotly(qualPlot, hoverinfo = c("Base", "Mean", "Median",
                                                           "Upper_Quartile", "Lower_Quartile",
                                                           "`10th_Percentile`","`90th_Percentile`"))
                )
                # Set the hoverinfo for bg rectangles to the vertices only,
                # This will effectively hide them
                qualPlot$x$data[[1]]$hoveron <- "points"
                qualPlot$x$data[[2]]$hoveron <- "points"
                qualPlot$x$data[[3]]$hoveron <- "points"
              }
            }

            if (plotType == "heatmap"){

              stopifnot(plotValue %in% c("Mean", "Median"))
              stopifnot(is.logical(cluster))

              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              if ("size" %in% names(dotArgs)){
                sz <- dotArgs$size
              }
              else{
                sz <- 0.2
              }
              if ("colour" %in% names(dotArgs) || "color" %in% names(dotArgs)){
                i <- which(names(dotArgs) %in% c("colour", "color"))
                lineCol <- dotArgs[[i]]
              }
              else{
                lineCol <- "grey20"
              }
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

              # Sort out the start positions
              df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))

              # Select the Mean or Median
              df <- df[c("Filename", "Start", plotValue, "Base")]

              #split data into correct lengths and fill NA's
              df <- split(df, f = df$Filename) %>%
                lapply(function(x){
                  Longest_sequence <- max(as.integer(gsub(".*-([0-9]*)", "\\1", x$Base)))
                  dfFill <- data.frame(Start = 1:Longest_sequence)
                  x <- dplyr::right_join(x, dfFill, by = "Start") %>%
                    zoo::na.locf()
                }) %>%
                dplyr::bind_rows()
              df$Start <- as.integer(df$Start)
              df <- df[!colnames(df) == c("Longest_sequence", "Base")]

              df <- reshape2::dcast(df, Filename ~ Start, value.var = plotValue)

              #cluster names true hclust names
              if(cluster){
                xx <- dplyr::select(df, -Filename)
                xx[is.na(xx)] <- 0
                clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
                row.ord <- order.dendrogram(clus)
                df <- df[row.ord,]
              }

              #melt data back to long format and change filenames to fit labels
              key <- df$Filename
              df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Start", value.name = "Data")
              df$Filename <- labels[df$Filename]

              # Reorganise the data frame
              df[[plotValue]] <- as.numeric(df$Data)
              df$Start <- as.integer(as.character(df$Start))
              df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))
              maxVal <- max(df[[plotValue]], na.rm = TRUE)
              phredMax <- ifelse(maxVal <= warn, 41, ceiling(maxVal + 1))

              # Start the heatmap
              qualPlot <- ggplot(df, aes_string(x = "Start", y = "Filename", fill = plotValue)) +
                geom_tile() +
                ngsReports:::scale_fill_pwf(na.omit(df[[plotValue]]), pwfCols, c(0, fail, warn, phredMax), FALSE, na.value = "white") +
                theme(panel.grid.minor = element_blank(),
                      panel.background = element_blank()) +
                scale_x_continuous(expand = c(0,0))

              if (usePlotly){

                qualPlot <- qualPlot +
                  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
                if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme

                t <- getSummary(x)
                t <- t[t$Category == "Per base sequence quality",]
                t$Filename <- labels[t$Filename]
                t$Filename <- factor(t$Filename, levels = levels(df$Filename))
                t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")

                sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)

                #plot dendrogram
                if(dendrogram && cluster){

                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0,0.5))

                  dendro <- suppressMessages(
                    plotly::ggplotly(dendro)
                  )

                  qualPlot <- suppressMessages(
                    plotly::subplot(dendro, sideBar, qualPlot, widths = c(0.1, 0.08, 0.82), margin = 0.001, shareY = TRUE) %>%
                      plotly::layout(xaxis3 = list(title = "Sequencing Cycle"))
                  )
                }
                else{
                  qualPlot <- suppressMessages(
                    plotly::subplot(plotly::plotly_empty(),
                                    sideBar, qualPlot,
                                    widths = c(0.1,0.08,0.82), margin = 0.001, shareY = TRUE) %>%
                      plotly::layout(xaxis3 = list(title = "Sequencing Cycle"),
                                     annotations = list(text = "Filename", showarrow = FALSE,
                                                        textangle = -90))
                    )
                }
              }
              else{
                # Add the custom themes
                if (!is.null(userTheme)) qualPlot <- qualPlot + userTheme
              }
            }

            qualPlot

          }
)
