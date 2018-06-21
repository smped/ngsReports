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
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
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
#' @importFrom scales rescale
#' @importFrom grDevices rgb
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @name plotSequenceContent
#' @rdname plotSequenceContent-methods
#' @export
setGeneric("plotSequenceContent",function(x, usePlotly = FALSE, labels, ...){standardGeneric("plotSequenceContent")})
#' @aliases plotSequenceContent,character
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "character",
          function(x, usePlotly = FALSE, labels, ...){
            x <- getFastqcData(x)
            plotSequenceContent(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceContent,FastqcFile
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, labels, ...){
            x <- getFastqcData(x)
            plotSequenceContent(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceContent,FastqcFileList
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, labels, ...){
            x <- getFastqcData(x)
            plotSequenceContent(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceContent,FastqcData
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, ...){

            # Get the SequenceContent
            df <- Per_base_sequence_content(x)
            
            if (!length(df)) {
              #stop("No sequence content Module")
              scPlot <- emptyPlot("No Sequence Content Module Detected")
              
              if(usePlotly) scPlot <- ggplotly(scPlot, tooltip = "")
              return(scPlot)
            }
            
            
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
            df <- reshape2::melt(df, id.vars = c("Filename", "Start"))
            colnames(df) <- c("Filename", "Position", "Base", "Percent")
            df$Base <- factor(df$Base, levels = c("T", "C", "A", "G"))
            df$Percent <- round(df$Percent, 2)

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
                plotly::subplot(plotly::plotly_empty(), scPlot, widths = c(0.14,0.86)) %>% 
                  layout(xaxis2 = list(title = "Position in read (bp)"), yaxis2 = list(title = "Percent (%)")))
            }

            scPlot

          }
)
#' @aliases plotSequenceContent,FastqcDataList
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, plotType = c("heatmap", "line"), pwfCols,
                   cluster = FALSE, dendrogram = TRUE, ...){

            # Get the SequenceContent
            df <- Per_base_sequence_content(x)
            
            if (!length(df)) {
              #stop("No sequence content Module")
              scPlot <- emptyPlot("No Sequence Content Module Detected")
              
              if(usePlotly) scPlot <- ggplotly(scPlot, tooltip = "")
              return(scPlot)
            }
            
            df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))
            df$End <- as.integer(gsub("[0-9]*-([0-9]*)", "\\1", df$Base))

            plotType <- match.arg(plotType)
            if (missing(pwfCols)) pwfCols <- ngsReports::pwf

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

            if (plotType == "heatmap"){

              df$A <- round(df$A, 2)
              df$C <- round(df$C, 2)
              df$G <- round(df$G, 2)
              df$T <- round(df$T, 2)
              maxBase <- max(vapply(c("A", "C", "G", "T"), function(x){max(df[[x]])}, numeric(1)))
              df$opacity <- 1 - df$G / maxBase
              df$colour <- with(df, rgb(red = `T` * opacity / maxBase,
                                        green = A * opacity / maxBase,
                                        blue = C * opacity / maxBase))

              basicStat <- Basic_Statistics(x)[c("Filename", "Longest_sequence")]

              df <- dplyr::right_join(df, basicStat, by = "Filename")
              df <- df[c("Filename", "Start", "End", "colour", "Longest_sequence", "A", "C", "G", "T")]

              # df <- split(df, f = df$Filename) %>%
              #   lapply(function(x){
              #     dfFill <- data.frame(Start = 1:x[["Longest_sequence"]][1])
              #     x <- dplyr::right_join(x, dfFill, by = "Start") %>%
              #       zoo::na.locf()
              #   }) %>%
              #   dplyr::bind_rows()
              df$Position <- as.integer(df$Start)
              df$Filename <- labels[df$Filename]

              if(cluster){
                mat <- reshape2::melt(df, id.vars = c("Filename", "Start"), measure.vars = c("A", "C", "G", "T"), variable.name = "Base", value.name = "Percent")
                
                mat <- reshape2::acast(mat, Filename ~ Start + Base, value.var = "Percent")
                mat[is.na(mat)] <- 0
                clus <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
                row.ord <- order.dendrogram(clus)
                df$Filename <- factor(df$Filename, levels = unique(df$Filename)[row.ord])
              }
              else{
                df$Filename <- factor(df$Filename, levels = unique(df$Filename))
              }
              
              key <-levels(df$Filename) 
              key <- names(labels[match(key, labels)])
              tileCols <- unique(df$colour)
              names(tileCols) <- unique(df$colour)
              
              df$ymax <- as.integer(df$Filename) + 0.5
              df$ymin <- df$ymax - 1
              df$xmax <- df$End + 0.5
              df$xmin <- df$Start - 1
              df$Window <- paste(df$Start, "-", df$End, "bp")
              
              scPlot <- ggplot(df, aes_string(fill = "colour", A = "A", C = "C", G = "G", T = "T", 
                                     Filename = "Filename", Window = "Window")) + 
                geom_rect(aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax")) +
                scale_fill_manual(values = tileCols) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0)) +
                theme_bw() +
                theme(legend.position = "none",
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank()) +
                labs(x = "Position in read (bp)",
                     y = "Filename")
              
              
              
              # scPlot <- ggplot(df, aes_string(x = "Position", y = "Filename", fill = "colour",
              #                                 A = "A", C = "C", G = "G", T = "T")) +
              #   geom_tile() +
              #   scale_fill_manual(values = tileCols) +
              #   scale_x_continuous(expand = c(0, 0)) +
              #   scale_y_discrete(expand = c(0, 0)) +
              #   theme_bw() +
              #   theme(legend.position = "none",
              #         panel.grid.minor = element_blank(),
              #         panel.grid.major = element_blank()) +
              #   labs(x = "Position in read (bp)",
              #        y = "Filename")
              
              if (!is.null(userTheme)) scPlot <- scPlot + userTheme
              
              if (usePlotly){
                scPlot <- scPlot +
                  theme(axis.ticks.y = element_blank(),
                        axis.text.y = element_blank())
                
                t <- getSummary(x)
                t <- t[t$Category == "Per base sequence content",]
                t$Filename <- labels[t$Filename]
                t$Filename <- factor(t$Filename, levels = levels(df$Filename))
                sideBar <- ngsReports:::makeSidebar(status = t, key = key, pwfCols = pwfCols)
                
                #plot dendro
                if (cluster && dendrogram){
                  dx <- ggdendro::dendro_data(clus)
                  dendro <- ngsReports:::ggdend(dx$segments) +
                    coord_flip() +
                    scale_y_reverse(expand = c(0, 0)) +
                    scale_x_continuous(expand = c(0, 0.5)) +
                    theme(panel.background = element_blank(),
                          panel.grid = element_blank())

                  
                  
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
                
                ## manually edit tooltip to remove colour
                sc <- lapply(1:length(scPlot$x$data), function(x){
                  scPlot$x$data[[x]]$text <<- gsub("colour:.*<br />A", "A",  scPlot$x$data[[x]]$text)
                  
                })
                
              }
            }
            if (plotType == "line"){
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
              # else{
              #   scPlot
              # }
              # scPlot
            }
           
            scPlot
          }
)
