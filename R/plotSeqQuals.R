#' @title Plot the Per Sequence Quality Scores
#'
#' @description Plot the Per Sequence Quality Scores for a set of FASTQC reports
#'
#' @details Plots the distribution of average sequence quality scores across the
#' set of files. Values can be plotted either as counts (`counts = TRUE`)
#' or as frequencies (`counts = FALSE`).
#'
#' Any faceting or scale adjustment can be performed after generation of the
#' initial plot, using the standard methods of ggplot2 as desired.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or path
#' @param counts `logical`. Plot the counts from each file if
#' `counts = TRUE`, otherwise the frequencies will be plotted
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All file names must be present in the names of the vector.
#' @param pattern Regex to remove from the end of any filenames
#' @param plotType `character`. Can only take the values
#' `plotType = "heatmap"` or `plotType = "line"`
#' @param warn,fail The default values for warn and fail are 5 and 10
#' respectively (i.e. percentages)
#' @param showPwf logical(1) Show PASS/WARN/FAIL status
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param colour Colour for single line plots
#' @param alpha set alpha for line graph bounds
#' @param heatCols Colour palette for the heatmap
#' @param heat_w Relative width of any heatmap plot components
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#' @param plotlyLegend logical(1) Show legend for interactive line plots
#' @param scaleFill,scaleColour ggplot2 scales
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # The default plot
#' plotSeqQuals(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fqName(fdl))
#' plotSeqQuals(fdl[r1])
#'
#' @docType methods
#'
#' @importFrom stats hclust dist
#' @importFrom scales percent comma percent_format
#' @importFrom tidyr pivot_wider complete nesting
#' @importFrom grDevices hcl.colors
#' @importFrom rlang !! sym
#' @import ggplot2
#'
#' @name plotSeqQuals
#' @rdname plotSeqQuals-methods
#' @export
setGeneric(
    "plotSeqQuals",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols, ...
    ){
        standardGeneric("plotSeqQuals")
    }
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod(
    "plotSeqQuals", signature = "ANY",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols, ...
    ){
        .errNotImp(x)
    }
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod(
    "plotSeqQuals", signature = "character",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols, ...
    ){
        x <- FastqcDataList(x)
        if (length(x) == 1) x <- x[[1]]
        plotSeqQuals(x, usePlotly, labels, pattern, pwfCols, ...)
    }
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod(
    "plotSeqQuals", signature = "FastqcData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
        showPwf = TRUE, counts = FALSE, alpha = 0.1, warn = 30, fail = 20,
        colour = "red", plotlyLegend = FALSE, ...
    ){

        df <- getModule(x, "Per_sequence_quality_scores")

        if (!length(df)) {
            p <- .emptyPlot("No Sequence Quality Moudule Detected")
            if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
            return(p)
        }

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern)
        labels <- labels[names(labels) %in% df$Filename]
        df$Filename <- labels[df$Filename]

        ## Sort out the colours
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))
        pwfCols <- setAlpha(pwfCols, alpha)
        stopifnot(warn > fail)

        ## Find the minimum quality value
        minQ <- min(df$Quality)

        ## make Ranges for rectangles and set alpha
        rects <- tibble(
            ymin = 0, ymax = max(df$Count),
            xmin = c(0, fail, warn), xmax = c(fail, warn, 41),
            Status = c("FAIL", "WARN", "PASS")
        )
        xLab <- "Mean Sequence Quality Per Read (Phred Score)"
        yLab <- "Number of Sequences"

        if (!counts) {

            ## Summarise to frequencies & initialise the plot
            df <- dplyr::group_by(df, Filename)
            df <- dplyr::mutate(df, Frequency = Count / sum(Count))
            df <- dplyr::ungroup(df)
            df$Frequency <- round(df$Frequency, 3)
            df$Percent <- scales::percent(df$Frequency)
            rects$ymax <- max(df$Frequency)

            p <- ggplot(df)
            if (showPwf) p <- p + geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                )
            )
            p <- p + geom_line(aes(Quality, Frequency), colour = colour)
            yLabelFun <- scales::percent_format(accuracy = 1)
            yLab <- "Frequency"

        }
        else{
            ## Initialise the plot using counts
            p <- ggplot(df)
            if (showPwf) p <- p + geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                )
            )
            p <- p + geom_line(aes(x = Quality, y = Count), colour = colour)
            yLabelFun <- scales::comma

        }

        p <- p +
            scale_fill_manual(values = getColours(pwfCols))  +
            scale_y_continuous(
                limits = c(0, rects$ymax[1]), expand = c(0, 0),
                labels = yLabelFun
            ) +
            scale_x_continuous(expand = c(0, 0)) +
            facet_wrap(~Filename) +
            labs(x = xLab, y = yLab) +
            theme_bw()
        p <- .updateThemeFromDots(p, ...)

        if (usePlotly) {

            if (!plotlyLegend)  p <- p + theme(legend.position = "none")
            ## Render as a plotly object
            hv <-  c("x", "y", "colour")
            p <- suppressMessages(
                suppressWarnings(plotly::ggplotly(p, hoverinfo = hv))
            )

            ## Set the hoverinfo for bg rectangles to the vertices only,
            ## This will effectively hide them
            p$x$data <- lapply(p$x$data, function(x){
                ## If there is a name component & it contains
                ## PASS/WARN/FAIL set the hoverinfo to none
                if ("name" %in% names(x)) {
                    if (grepl("(PASS|WARN|FAIL)", x$name)) {
                        x$hoverinfo <- "none"
                    }
                }
                x
            })
        }

        ## Draw the plot
        p

    }
)
#' @rdname plotSeqQuals-methods
#' @export
setMethod(
    "plotSeqQuals", signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
        counts = FALSE, alpha = 0.1, warn = 30, fail = 20, showPwf = TRUE,
        plotType = c("heatmap", "line"), dendrogram = FALSE, cluster = FALSE,
        scaleFill = NULL, heatCols = hcl.colors(100, "inferno"), heat_w = 8,
        scaleColour = NULL, plotlyLegend = FALSE, ...
    ){

        ## Read in data
        mod <- "Per_sequence_quality_scores"
        df <- getModule(x, "Per_sequence_quality_scores")

        if (!length(df)) {
            p <- .emptyPlot("No Sequence Quality Moudule Detected")
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        ## As the Quality Scores may have different ranges across different
        ## samples, fill missing values with zero
        df <- complete(df, Filename, nesting(Quality), fill = list(Count = 0))
        df$Quality <- as.integer(df$Quality)

        ## Check for valid plotType
        plotType <- match.arg(plotType)
        xLab <- "Mean Sequence Quality Per Read (Phred Score)"
        plotVal <- ifelse(counts, "Count", "Frequency")
        stopifnot(warn > fail)

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern)
        labels <- labels[names(labels) %in% df$Filename]

        ## Sort out the colours
        if (base::missing(pwfCols)) pwfCols <- pwf

        if (plotType == "heatmap") {

            ## Summarise to frequencies & initialise the plot
            df <- dplyr::group_by(df, Filename)
            df <- dplyr::mutate(df, Frequency = Count / sum(Count))
            df <- dplyr::ungroup(df)
            df$Frequency <- round(df$Frequency, 3)
            df$Total <- scales::comma(df$Count)

            ## Now define the order for a dendrogram if required
            ## This only applies to a heatmap
            key <- names(labels)
            cols <- c("Filename", "Quality", "Frequency")
            clusterDend <- .makeDendro(
                df[cols], "Filename","Quality", "Frequency"
            )
            dx <- ggdendro::dendro_data(clusterDend)
            if (dendrogram | cluster) key <- labels(clusterDend)
            df$Filename <- factor(labels[df$Filename], levels = labels[key])
            if (!dendrogram) dx$segments <- dx$segments[0,]

            if (is.null(scaleFill)) {
                scaleFill <- scale_fill_viridis_c(
                    option = "inferno", limits = c(0, 1)
                )
                if (!is.null(heatCols)) {
                    scaleFill <- scale_fill_gradientn(
                        colours = heatCols, limits = c(0, 1)
                    )
                }
            }
            stopifnot(is(scaleFill, "ScaleContinuous"))
            stopifnot(scaleFill$aesthetics == "fill")

            hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
            p <- ggplot(df, aes(Quality, Filename, C = Total)) +
                geom_tile(aes(fill = Frequency)) +
                ggtitle("Per Sequence Quality Score") +
                labs(x = xLab, y = c()) +
                scaleFill +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0), position = "right") +
                theme(
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    plot.title = element_text(hjust = hj),
                    axis.title.x = element_text(hjust = hj)
                )
            if (showPwf | dendrogram)
                p <- p +
                theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))
            p <- .updateThemeFromDots(p, ...)

            status <- getSummary(x)
            status <- subset(status, Category == gsub("_", " ", mod))
            status$Filename <- factor(
                labels[status$Filename], levels = labels[key]
            )
            if (!showPwf) status <- status[0, ]

            p <- .prepHeatmap(
                p, status, dx$segments, usePlotly, heat_w, pwfCols
            )

        }

        if (plotType == "line") {

            ## make Ranges for rectangles and set alpha
            pwfCols <- setAlpha(pwfCols, alpha)
            rects <- tibble(
                ymin = 0, ymax = max(df$Count),
                xmin = c(0, fail, warn), xmax = c(fail, warn, 41),
                Status = c("FAIL", "WARN", "PASS")
            )
            ## No clustering required so just use the labels
            df$Filename <- labels[df$Filename]

            yLab <- ifelse(counts, "Number of Sequences", "Frequency")
            yLabelFun <- ifelse(
                counts, scales::comma, scales::percent_format(accuracy = 1)
            )
            plotVal <- ifelse(counts, "Count", "Frequency")

            if (!counts) {
                ## Summarise to frequencies & initialise the plot
                df <- dplyr::group_by(df, Filename)
                df <- dplyr::mutate(df, Frequency = Count / sum(Count))
                df <- dplyr::ungroup(df)
                df$Frequency <- round(df$Frequency, 4)
                rects$ymax <- max(df$Frequency)
            }

            if (!is.null(scaleColour)) {
                stopifnot(is(scaleColour, "ScaleDiscrete"))
                stopifnot(scaleColour$aesthetics == "colour")
            }

            p <- ggplot(df)
            if (showPwf) p <- p + geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                )
            ) +
                scale_fill_manual(values = getColours(pwfCols))

            p <- p +
                geom_line(
                    aes(x = Quality, y = !!sym(plotVal), colour = Filename)
                ) +
                scale_y_continuous(
                    limits = c(0, rects$ymax[1]), expand = c(0, 0),
                    labels = yLabelFun
                ) +
                scale_x_continuous(expand = c(0, 0)) +
                scaleColour +
                labs(x = xLab, y = yLab) +
                guides(fill = "none") +
                theme_bw()
            p <- .updateThemeFromDots(p, ...)

            if (usePlotly) {

                hv <- c("x", "y", "colour")
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                p <- suppressMessages(
                    suppressWarnings(plotly::ggplotly(p, hoverinfo = hv))
                )

                ## Turn off the hoverinfo for the bg rectangles
                ## This will effectively hide them
                p$x$data <- lapply(p$x$data, function(x){
                    ## If there is a name component & it contains
                    ## PASS/WARN/FAIL set the hoverinfo to none
                    if ("name" %in% names(x)) {
                        if (grepl("(PASS|WARN|FAIL)", x$name))
                            x$hoverinfo <- "none"
                    }
                    x
                })
            }}

        p
    }

)
