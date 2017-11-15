#' @title Plot the per base content as a heatmap
#'
#' @description Plot the Per Base content for a set of FASTQC files.
#' Informative plot where per base sequence content (%A, %T, %G, %C),
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly \code{logical}. Generate an interactive plot using plotly
#' @param plotType \code{character}. Type of plot to generate. Must be "line" or "heatmap"
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A ggplot2 object
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
#' plotSequenceContent(fdl)
#'
#'
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom scales rescale
#' @importFrom grDevices rgb
#' @importFrom magrittr %>%
#'
#' @name plotSequenceContent
#' @rdname plotSequenceContent-methods
#' @export
setGeneric("plotSequenceContent",function(x, usePlotly = FALSE, ...){standardGeneric("plotSequenceContent")})
#' @aliases plotSequenceContent,character
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceContent(x, usePlotly,...)
          }
)
#' @aliases plotSequenceContent,FastqcFile
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceContent(x, usePlotly,...)
          }
)
#' @aliases plotSequenceContent,FastqcFileList
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotSequenceContent(x, usePlotly,...)
          }
)
#' @aliases plotSequenceContent,FastqcData
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, ...){

            # Get the SequenceContent
            df <- Per_base_sequence_content(x)
            df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))

            # Drop the suffix, or check the alternate labels
            if (missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
            }
            else{
              if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            df$Filename <- labels[df$Filename]
            df <- df[!colnames(df) == "Base"]
            df <- melt(df, id.vars = c("Filename", "Start"))
            colnames(df) <- c("Filename", "Position", "Base", "Percent")
            df$Base <- factor(df$Base, levels = c("T", "C", "A", "G"))

            #set colours
            baseCols <- c(`T`="red", G = "black", A = "green", C = "blue")

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            scPlot <- ggplot(df, aes_string(x = "Position", y = "Percent", colour = "Base")) +
              geom_line() +
              facet_wrap(~Filename) +
              scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
              scale_x_continuous(expand = c(0, 0)) +
              guides(fill = FALSE) +
              labs(x = "Position in read (bp)",
                   y = "Percent (%)") +
              theme_bw() +
              scale_colour_manual(values = baseCols)
            if(usePlotly){

              scPlot <- suppressMessages(
                subplot(plotly::plotly_empty(), scPlot, widths = c(0.08,0.92), margin = 0.001) %>%
                  layout(xaxis2 = list(title = "Position in read (bp)"))
              )
            }

            scPlot

          }
)
#' @aliases plotSequenceContent,FastqcDataList
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, plotType = "heatmap", pwfCols,
                   clusterNames = FALSE, dendrogram = TRUE, ...){

            stopifnot(plotType %in% c("heatmap", "line"))
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf

            # Get the SequenceContent
            df <- Per_base_sequence_content(x)
            df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))

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

            if(plotType == "heatmap"){

              maxBase <- max(vapply(c("A", "C", "G", "T"), function(x){max(df[[x]])}, numeric(1)))
              df$opacity <- 1 - df$G / maxBase
              df$colour <- with(df, rgb(red = round(`T` * opacity / maxBase, 2),
                                        green = round(A * opacity / maxBase, 2),
                                        blue = round(C * opacity / maxBase, 2)))

              basicStat <- Basic_Statistics(x)[c("Filename", "Longest_sequence")]

              df <- dplyr::right_join(df, basicStat, by = "Filename")
              df <- df[c("Filename", "Start", "colour", "Longest_sequence", "A", "C", "G", "T")]

              df <- split(df, f = df$Filename) %>%
                lapply(function(x){
                  dfFill <- data.frame(Start = 1:x[["Longest_sequence"]][1])
                  x <- dplyr::right_join(x, dfFill, by = "Start") %>%
                    zoo::na.locf()
                }) %>%
                dplyr::bind_rows()
              df$Position <- as.integer(df$Start)
              df$Filename <- labels[df$Filename]

              if(clusterNames){
                mat <- reshape2::melt(df, id.vars = c("Filename", "Start"), measure.vars = c("A", "C", "G", "T"), variable.name = "Base", value.name = "Percent")
                mat <- reshape2::acast(mat, Filename ~ Start + Base, value.var = "Percent")
                clus <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
                row.ord <- order.dendrogram(clus)
                df$Filename <- factor(df$Filename, levels = unique(df$Filename)[row.ord])
              }
              else{
                df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))
              }

              key <- levels(df$Filename)
              tileCols <- unique(df$colour)
              names(tileCols) <- unique(df$colour)

              scPlot <- ggplot(df, aes_string(x = "Position", y = "Filename", fill = "colour",
                                              A = "A", C = "C", G = "G", T = "T")) +
                geom_tile() +
                scale_fill_manual(values = tileCols) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0)) +
                theme(legend.position = "none",
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank()) +
                labs(x = "Position in read (bp)",
                     y = "Filename")

              if (!is.null(userTheme)) scPlot <- scPlot + userTheme

              if (usePlotly){
                scPlot <- scPlot +
                  theme(axis.ticks.y = element_blank(),
                        axis.text.y = element_blank())

                t <- getSummary(x)
                t <- t[t$Category == "Per base sequence content",]
                t$Filename <- labels[t$Filename]
                t$Filename <- factor(t$Filename, levels = levels(df$Filename))
                sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)

                #plot dendro
                if (clusterNames && dendrogram){
                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0, 0.5)) +
                    theme(panel.background = element_rect(fill = "white"))

                  scPlot <- suppressWarnings(
                    suppressMessages(
                      plotly::subplot(dendro, sideBar, scPlot, widths = c(0.1,0.08,0.82),
                                      margin = 0.001, shareY = TRUE)
                    ))
                }
                else{
                  # Return the plot
                  scPlot <- suppressWarnings(
                    suppressMessages(
                      plotly::subplot(plotly::plotly_empty(), sideBar, scPlot, widths = c(0.1,0.08,0.82),
                                      margin = 0.001, shareY = TRUE)
                    )
                  )
                }

              }
              else{
                scPlot
              }
              scPlot
            }
            else{
              df$Filename <- labels[df$Filename]
              df <- df[!colnames(df) == "Base"]
              df <- melt(df, id.vars = c("Filename", "Start"))
              colnames(df) <- c("Filename", "Position", "Base", "Percent")
              df$Base <- factor(df$Base, levels = c("T", "C", "A", "G"))

              #set colours
              baseCols <- c(`T`="red", G = "black", A = "green", C = "blue")

              scPlot <- ggplot(df, aes_string(x = "Position", y = "Percent", colour = "Base")) +
                geom_line() +
                facet_wrap(~Filename) +
                scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                scale_x_continuous(expand = c(0, 0)) +
                guides(fill = FALSE) +
                labs(x = "Position in read (bp)",
                     y = "Percent (%)") +
                theme_bw() +
                scale_colour_manual(values = baseCols)

              if (!is.null(userTheme)) scPlot <- scPlot + userTheme

              if(usePlotly){
                scPlot <- suppressMessages(
                  ggplotly(scPlot) %>% layout(legend = list(x = 0.85, y = 1))
                )
              }
              else{
                scPlot
              }
              scPlot
            }
          }
)
