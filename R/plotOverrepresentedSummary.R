#' @title Plot a summary of Over-represented Sequences
#'
#' @description Plot a summary of Over-represented Sequences for a set of FASTQC reports
#'
#' @details Percentages are obtained by simply summing those within a report.
#' Any possible double counting by FastQC is ignored for the purposes of a simple approximation.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param n The number of sequences to plot from an individual file
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param paletteName Name of the palette for colouring the possible sources of the overrepresented sequences.
#' Must be a palette name from \code{RColorBrewer}.
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
#' # Another example which isn't ideal
#' plotOverrepresentedSummary(fdl)
#'
#' @importFrom tidyr spread
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar geom_text
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlim ylim
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme_bw theme_void
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 facet_grid
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom plotly layout
#' @importFrom grDevices rgb
#'
#'
#' @name plotOverrepresentedSummary
#' @rdname plotOverrepresentedSummary-methods
#' @export
setGeneric("plotOverrepresentedSummary",function(x, usePlotly = FALSE, ...){standardGeneric("plotOverrepresentedSummary")})
#' @aliases plotOverrepresentedSummary,character
#' @rdname plotOverrepresentedSummary-methods
#' @export
setMethod("plotOverrepresentedSummary", signature = "character",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotOverrepresentedSummary(x, usePlotly,...)
          }
)
#' @aliases plotOverrepresentedSummary,FastqcFile
#' @rdname plotOverrepresentedSummary-methods
#' @export
setMethod("plotOverrepresentedSummary", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotOverrepresentedSummary(x, usePlotly,...)
          }
)
#' @aliases plotOverrepresentedSummary,FastqcFileList
#' @rdname plotOverrepresentedSummary-methods
#' @export
setMethod("plotOverrepresentedSummary", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
            x <- getFastqcData(x)
            plotOverrepresentedSummary(x, usePlotly,...)
          }
)
#' @aliases plotOverrepresentedSummary,FastqcData
#' @rdname plotOverrepresentedSummary-methods
#' @export
setMethod("plotOverrepresentedSummary", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, n = 10, pwfCols, ...){

            df <- Overrepresented_sequences(x)

            if (nrow(df) == 0) {
              #stop("No overrepresented sequences were detected by FastQC")
              overPlot <- ggplot() +
                geom_text(aes(x = 0.5, y = 0.8, label = "No overrepresented sequences")) +
                theme_void() +
                xlim(c(0, 1)) +
                ylim(c(0, 1))
              return(overPlot)
            }

            # Drop the suffix, or check the alternate labels
            if(base::missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
            df$Filename <- labels[df$Filename]

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            df <- dplyr::top_n(df, n, Percentage)
            df$Status <- cut(df$Percentage, breaks = c(0, 0.1, 1, 100), labels = c("PASS", "WARN", "FAIL"))
            df$Possible_Source <- gsub(" \\([0-9]*\\% over [0-9]*bp\\)", "", df$Possible_Source)
            df$Sequence <- factor(df$Sequence, levels = rev(df$Sequence))
            df <- droplevels(df)
            ymax <- 1.05*max(df$Percentage)

            # Sort out the colours & pass/warn/fail breaks
            if (missing(pwfCols)) pwfCols <- getColours(ngsReports::pwf)
            pwfCols <- pwfCols[names(pwfCols) %in% levels(df$Status)]

            if (usePlotly){
              overPlot <- plot_ly(df, x = ~Percentage, y = ~Sequence, type = "bar",
                                  name = ~Possible_Source, color = I(pwfCols[as.character(df$Status)]),
                                  hoverinfo = "text",
                                  text = ~paste("Filename: ",
                                                Filename,
                                                "<br> Percentage: ",
                                                Percentage,
                                                "<br> Sequence: ",
                                                Sequence,
                                                "<br> Possible Source: ",
                                                Possible_Source)) %>%
                layout(xaxis = list(title = "Percent of Total Reads"),
                       yaxis = list(title = "Overrepresented Sequence",
                                    showticklabels = FALSE), title = df$Filename[1])
            }
            else{
              overPlot <- ggplot(df, aes_string(x = "Sequence", y = "Percentage", fill = "Status")) +
                geom_bar(stat = "identity") +
                ylab("Overrepresented Sequences (% of Total)") +
                scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
                theme_bw() +
                coord_flip() +
                facet_grid(Possible_Source~., scales = "free_y", space = "free") +
                scale_fill_manual(values = pwfCols)

              # Add the basic customisations
              if (!is.null(userTheme)) overPlot <- overPlot + userTheme

            }

            overPlot

            # Add functionality to export a FASTA file of the sequences to the shiny app?
            # This will obviously work best under plotly as the names will be silly otherwise

          }
)
#' @aliases plotOverrepresentedSummary,FastqcDataList
#' @rdname plotOverrepresentedSummary-methods
#' @export
setMethod("plotOverrepresentedSummary", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, ..., paletteName = "Set1"){

            df <- Overrepresented_sequences(x)

            if (nrow(df) == 0) stop("No overrepresented sequences were detected by FastQC")

            # Drop the suffix, or check the alternate labels
            if(base::missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            # Get any arguments for dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(ggplot2::theme))
            keepArgs <- which(names(dotArgs) %in% allowed)
            userTheme <- c()
            if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

            Possible_Source <- c() # Here to avoid a NOTE in R CMD check...
            df$Possible_Source <- gsub(" \\([0-9]*\\% over [0-9]*bp\\)", "", df$Possible_Source)
            df <- dplyr::group_by(df, Filename, Possible_Source)
            df <- dplyr::summarise(df, Percentage = sum(Percentage))
            maxChar <- max(nchar(df$Filename))
            df$Filename <- labels[df$Filename]
            df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

            # Set the axis limits. Just scale the upper limit by 1.05
            ymax <- max(dplyr::summarise(dplyr::group_by(df, Filename),
                                         Total = sum(Percentage))$Total)*1.05

            # Define the palette
            stopifnot(paletteName %in% rownames(RColorBrewer::brewer.pal.info))
            nSource <- length(unique(df$Possible_Source))
            pal <- RColorBrewer::brewer.pal(nSource, paletteName) %>% set_names(unique(df$Possible_Source))

            if (usePlotly){

              #set left margin
              if(maxChar < 10) l <- 80
              if(maxChar >= 10 & maxChar < 15) l <- 110
              if(maxChar >= 15 & maxChar < 20) l <- 130
              if(maxChar >= 20)  l <- 150


              overPlot <- plot_ly(df, x = ~Percentage, y = ~Filename, type = "bar", color = ~Possible_Source,
                      colors = pal, hoverinfo = "text",
                      text = ~paste("Filename: ",
                                    Filename,
                                    "<br> Percentage: ",
                                    Percentage,
                                    "<br> Possible Source: ",
                                    Possible_Source)) %>%
                layout(xaxis = list(title = "Percent of Total Reads"),
                       margin=list(l=l), barmode = "stack", showlegend = FALSE)

            }
            else{
              overPlot <- ggplot(df, aes_string(x = "Filename", y = "Percentage", fill = "Possible_Source")) +
                geom_bar(stat = "identity") +
                labs(y = "Overrepresented Sequences (% of Total)",
                     fill = "Possible Source") +
                scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
                scale_fill_manual(values = pal) +
                theme_bw() +
                coord_flip()

              # Add the basic customisations
              if (!is.null(userTheme)) overPlot <- overPlot + userTheme
            }

            overPlot

          }
)
