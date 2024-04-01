#' @title Plot the combined Sequence_Duplication_Levels information
#'
#' @description Plot the Sequence_Duplication_Levels information for a set of
#' FASTQC reports
#'
#' @details
#' This extracts the Sequence_Duplication_Levels from the supplied object and
#' generates a ggplot2 object, with a set of minimal defaults. For multiple
#' reports, this defaults to a heatmap with block sizes proportional to the
#' percentage of reads belonging to that duplication category.
#'
#' If setting `usePlotly = FALSE`, the output of this function can be
#' further modified using standard ggplot2 syntax. If setting
#' `usePlotly = TRUE` an interactive plotly object will be produced.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file path
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail values in the plot
#' @param warn,fail The default values for warn and fail are 20 and 50
#' respectively (i.e. percentages)
#' @param lineCol,lineWidth Colours and width of lines drawn
#' @param deduplication Plot Duplication levels 'pre' or 'post' deduplication.
#' Can only take values "pre" and "post"
#' @param plotType Choose between "heatmap" and "line"
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` Plot will automatically be clustered if TRUE.
#' @param heatCol Colour palette used for the heatmap
#' @param barFill,barCol Colours for bars when calling geom_col()
#' @param pattern regex to remove from the end of fastp & fastq file names
#' @param showPwf logical(1) Show PWF rectangles in the background
#' @param plotlyLegend logical(1) Show legend for line plots when using
#' interactive plots
#' @param scaleFill Discrete scale used to fill heatmap cells
#' @param plotTheme \link[ggplot2]{theme} object. Applied after a call to
#' theme_bw()
#' @param maxLevel The maximum duplication level to plot. Beyond this level, all
#' values will be summed
#' @param heat_w Relative width of the heatmap relative to other plot components
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#'
#' @return A standard ggplot2 or plotly object
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
#' # Draw the default plot for a single file
#' plotDupLevels(fdl[[1]])
#'
#' plotDupLevels(fdl)
#'
#' @docType methods
#'
#' @import ggplot2
#' @importFrom stats hclust dist
#' @importFrom grDevices hcl.colors
#' @importFrom tidyselect contains
#' @importFrom forcats fct_inorder
#' @import tibble
#' @importFrom rlang !! sym
#'
#' @name plotDupLevels
#' @rdname plotDupLevels-methods
#' @export
setGeneric(
    "plotDupLevels",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...) standardGeneric("plotDupLevels")
)
#' @rdname plotDupLevels-methods
#' @export
setMethod(
    "plotDupLevels", signature = "ANY",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){
        .errNotImp(x)
    }
)
#' @rdname plotDupLevels-methods
#' @export
setMethod(
    "plotDupLevels", signature = "FastqcData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        pwfCols, warn = 20, fail = 50, showPwf = TRUE, plotlyLegend = FALSE,
        lineCol = c("red", "blue"), lineWidth = 1, ...
    ){

        mod <- "Sequence_Duplication_Levels"
        df <- getModule(x, mod)

        if (!length(df)) {
            p <- .emptyPlot("No Duplication Levels Module Detected")
            if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
            return(p)
        }

        ## Convert from wide to long
        df <- tidyr::pivot_longer(
            df, cols = contains("Percentage"),
            names_to = "Type", values_to = "Percentage"
        )
        df$Duplication_Level <- forcats::fct_inorder(df$Duplication_Level)
        df$Type <- stringr::str_replace(df$Type, "Percentage_of_", "% ")
        df$Type <- stringr::str_to_title(df$Type)
        df$Type <- paste(df$Type, "sequences")
        df$x <- as.integer(df$Duplication_Level)
        df$Percentage <- round(df$Percentage, 2)

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern)
        labels <- labels[names(labels) %in% df$Filename]
        df$Filename <- labels[df$Filename]

        ## Sort out the colours
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))
        pwfCols <- setAlpha(pwfCols, 0.1)

        ## Set the background rectangles
        rects <- tibble(
            xmin = 0.5, xmax = max(df$x) + 0.5,
            ymin = c(0, warn, fail), ymax = c(warn, fail, 100),
            Status = c("PASS", "WARN", "FAIL")
        )

        ##Axis labels
        xlab <- gsub("_", " ", mod)
        ylab <- "Percentage (%)"

        p <- ggplot(data = df)
        if (showPwf) p <- p +
            geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                )
            )

        p <- p + geom_line(
            aes(x, Percentage, colour = Type, group = Type),
            linewidth = lineWidth
        ) +
            scale_fill_manual(values = getColours(pwfCols)) +
            scale_colour_manual(values = lineCol) +
            scale_x_continuous(
                breaks = unique(df$x),
                labels = levels(df$Duplication_Level), expand = c(0, 0)
            ) +
            scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
            facet_wrap(~Filename) +
            labs(x = xlab, y = ylab, colour = c()) +
            guides(fill = "none") +
            theme_bw() +
            theme(
                legend.position = "inside", legend.position.inside = c(1, 1),
                legend.justification.inside = c(1, 1),
                legend.background = element_rect(
                    colour = "black", linewidth = 0.2
                )
            )
        p <- .updateThemeFromDots(p, ...)

        if (usePlotly) {
            tt <- c("colour", "Percentage")
            if (!plotlyLegend) p <- p + theme(legend.position = "none")
            p <- suppressMessages(plotly::ggplotly(p, tooltip = tt))
            ## Make sure there are no hovers over the background rectangles
            p$x$data <- lapply(
                p$x$data,
                function(x) {
                    if (all(x$text == "")) x$hoverinfo <- "none"
                    x
                }
            )
        }

        p

    }
)
#' @rdname plotDupLevels-methods
#' @export
setMethod(
    "plotDupLevels",signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        pwfCols, warn = 20, fail = 50, showPwf = TRUE, plotlyLegend = FALSE,
        deduplication = c("pre", "post"), plotType = c("heatmap", "line"),
        cluster = FALSE, dendrogram = FALSE,
        heatCol = hcl.colors(50, "inferno"), heat_w = 8, ...
    ){

        mod <- "Sequence_Duplication_Levels"
        df <- getModule(x, mod)

        if (!length(df)) {
            p <- .emptyPlot("No Duplication Levels Module Detected")
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        ## Select the 'pre/post' option & clean up the data
        deduplication <- match.arg(deduplication)
        type <- dplyr::case_when(
            deduplication == "pre" ~ "Percentage_of_total",
            deduplication == "post" ~ "Percentage_of_deduplicated"
        )
        df <- df[c("Filename", "Duplication_Level", type)]
        df[[type]] <- round(df[[type]], 2)
        Duplication_Level <- c() # Here to avoid R CMD check issues

        ## Check the plotType
        plotType <- match.arg(plotType)

        ## These will begin in order, but may not stay this way
        ## in the following code
        df$Duplication_Level <- forcats::fct_inorder(df$Duplication_Level)
        dupLevels <- levels(df$Duplication_Level)

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern)
        labels <- labels[names(labels) %in% df$Filename]
        key <- names(labels)

        ## Sort out the colours
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))

        if (plotType == "line") {

            ## Now set everything as factors
            df$Filename <- forcats::fct_inorder(labels[df$Filename])
            df$x <- as.integer(df$Duplication_Level)

            ## Make transparent for a line plot
            pwfCols <- setAlpha(pwfCols, 0.1)

            ## Set the background rectangles
            rects <- tibble(
                xmin = 0.5,
                xmax = max(df$x) + 0.5,
                ymin = c(0, warn, fail),
                ymax = c(warn, fail, 100),
                Status = c("PASS", "WARN", "FAIL")
            )

            p <- ggplot(df, aes(label = !!sym("Duplication_Level")))

            if (showPwf) p <- p + geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                )
            )

            p <- p +
                geom_line(aes(x, y = !!sym(type), colour = Filename)) +
                scale_fill_manual(values = getColours(pwfCols)) +
                scale_x_continuous(
                    breaks = seq_along(dupLevels),labels = dupLevels,
                    expand = c(0, 0)
                ) +
                scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                labs(x = "Duplication Level", y = "Percentage of Total") +
                guides(fill = "none") +
                theme_bw()
            p <- .updateThemeFromDots(p, ...)

            if (usePlotly) {

                tt <- c("colour", type, "Duplication_Level")
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                p <- suppressMessages(plotly::ggplotly(p, tooltip = tt))

                ## Make sure there are no hovers over the background rectangles
                p$x$data <- lapply(p$x$data, .hidePWFRects)

            }


        }

        if (plotType == "heatmap") {

            dl <- "Duplication_Level"
            pt <- "Percentage_of_total"
            ## Set up the dendrogram
            clusterDend <- .makeDendro(df, "Filename", dl, type)
            dx <- ggdendro::dendro_data(clusterDend)
            if (dendrogram | cluster) {
                key <- labels(clusterDend)
            }
            if (!dendrogram) dx$segments <- dx$segments[0,]
            df$Filename <- factor(labels[df$Filename], levels = labels[key])

            ## Setup to plot in tiles for easier plotly compatability
            df <-  dplyr::arrange(df, Filename, Duplication_Level)
            df <- split(df, f = df[["Filename"]])
            df <- lapply(df, function(x){
                x$xmax <- cumsum(x[[pt]])
                x$xmax <- round(x[["xmax"]], 1) # Deal with rounding errors
                x$xmin <- c(0, x[["xmax"]][-nrow(x)])
                x
            })
            df <-  dplyr::bind_rows(df)
            df$ymax <- as.integer(df[["Filename"]]) + 0.5
            df$ymin <- df[["ymax"]] - 1

            ## Setup some more plotting parameters
            cols <- colorRampPalette(heatCol)(length(dupLevels))
            hj <-  0.5 * heat_w / (heat_w + 1 + dendrogram)
            p <- ggplot(
                df,
                aes(fill = !!sym(dl), total = !!sym(pt), label = Filename)
            ) +
                geom_rect(
                    aes(
                        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        colour = !!sym(dl)
                    )
                ) +
                ggtitle(gsub("_", " ", mod)) +
                scale_fill_manual(values = cols) +
                scale_colour_manual(values = cols) +
                scale_y_continuous(
                    breaks = seq_along(levels(df$Filename)),
                    labels = levels(df$Filename),
                    expand = c(0, 0), position = "right"
                ) +
                scale_x_continuous(expand = c(0, 0)) +
                labs(x = pt, fill = "Duplication\nLevel") +
                guides(colour = "none") +
                theme_bw() +
                theme(
                    plot.title = element_text(hjust = hj),
                    axis.title.x = element_text(hjust = hj)
                )
            if (dendrogram | showPwf)
                p <- p +
                theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))
            p <- .updateThemeFromDots(p, ...)

            status <- getSummary(x)
            status <- subset(status, Category == "Sequence Duplication Levels")
            status$Filename <-
                factor(labels[status$Filename], levels = levels(df$Filename))
            if (!showPwf) status <- status[0, ]
            p <- .prepHeatmap(
                p, status, dx$segments, usePlotly, heat_w, pwfCols
            )

        }

        p
    }
)
#' @importFrom stats weighted.mean
#' @importFrom rlang sym !!
#' @importFrom dplyr group_by ungroup summarise
#' @importFrom tidyr unnest
#' @importFrom scales percent
#' @rdname plotDupLevels-methods
#' @export
setMethod(
    "plotDupLevels", signature = "FastpData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        pwfCols, warn = 20, fail = 50, showPwf = FALSE, maxLevel = 10,
        lineCol = "red", barFill = "dodgerblue4", barCol = barFill,
        plotlyLegend = FALSE, plotTheme = theme_get(), ...
    ){

        mod <- "Duplication"
        df <- getModule(x, mod)

        ## Extract the histogram data
        df <- unnest(df, !!sym("histogram"))
        if (!nrow(df)) {
            msg <- "No Duplication Histogram Data Found"
            message(msg)
            p <- .emptyPlot(msg)
            if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
            return(p)
        }

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(df, labels, pattern = pattern)
        labels <- labels[names(labels) %in% df$Filename]
        df$Filename <- labels[df$Filename]

        ## Prepare for plotting
        df <- dplyr::arrange(df, !!sym("duplication_level"))
        df$duplication_level <- ifelse(
            df$duplication_level > maxLevel,
            paste0(">", maxLevel),
            as.character(df$duplication_level)
        )
        df$duplication_level <- forcats::fct_inorder(df$duplication_level)
        df <- group_by(df, Filename, !!sym("duplication_level"))
        df <- summarise(
            df,
            rate = unique(!!sym("rate")),
            histogram = sum(!!sym("histogram")),
            mean_gc = weighted.mean(
                !!sym("mean_gc"), w = !!sym("duplication_rate")
            ),
            duplication_rate = sum(!!sym("duplication_rate")), .groups = "drop"
        )
        df$x <- as.integer(df$duplication_level)
        df[["Duplication Level"]] <- df$duplication_level
        df[["Duplication Rate"]] <- percent(df$duplication_rate, 0.1)
        df[["Mean GC"]] <- percent(df$mean_gc, 0.1)
        main <- sprintf(
            "%s: Duplication Rate (%s)",
            unique(df$Filename), percent(unique(df$rate), 0.1)
        )

        ## Sort out the colours for any rectangles
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))
        pwfCols <- setAlpha(pwfCols, 0.2)
        stopifnot(is(plotTheme, "theme"))

        ## Set the limits & rectangles
        ylim <- c(0, max(df$duplication_rate) * 1.05)
        expand_x <- round(0.015*(max(df$x) - min(df$x)), 1)
        rects <- tibble(
            xmin = min(df$x - 0.5) - expand_x,
            xmax = max(df$x + 0.5) + expand_x,
            ymin = c(0, warn, fail) / 100,
            ymax = c(warn, fail, max(ylim * 100)) / 100,
            Status = forcats::fct_inorder(c("PASS", "WARN", "FAIL"))
        )
        xlim <- c(rects$xmin[[1]], rects$xmax[[1]])
        rects <- dplyr::filter(rects, ymin < max(ylim))
        rects$ymax[rects$ymax > max(ylim)] <- max(ylim)

        p <- ggplot(
            df,
            aes(
                xlab = !!sym("Duplication Level"),
                rate = !!sym("Duplication Rate"), gc = !!sym("Mean GC")
            )
        )
        if (showPwf) {
            p <- p + geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                ),
                inherit.aes = FALSE
            )
        }

        p <- p +
            geom_col(
                aes(x, !!sym("duplication_rate")), fill = barFill,
                colour = barCol
            ) +
            geom_line(
                aes(x, y = !!sym("mean_gc"), group = 1), col = lineCol[[1]], ...
            ) +
            scale_x_continuous(
                breaks = df$x, labels = levels(df$duplication_level),
                expand = c(0, 0), limits = xlim
            ) +
            scale_y_continuous(
                labels = percent, expand = c(0, 0), limits = ylim
            ) +
            scale_fill_manual(values = getColours(pwfCols)) +
            labs(x = "Duplication Level", y = "Duplication Rate") +
            ggtitle(main) +
            theme_bw() +
            plotTheme

        if (usePlotly) {
            if (!plotlyLegend) p <- p + theme(legend.position = "none")
            tt <- c("Duplication Level", "Duplication Rate", "Mean GC")
            p <- plotly::ggplotly(p, tooltip = tt)
            p$x$data <- lapply(p$x$data, .hidePWFRects)
        }

        p
    }
)
#' @importFrom tidyr unnest
#' @importFrom rlang sym !!
#' @importFrom dplyr group_by mutate summarise ungroup
#' @importFrom grDevices colorRampPalette
#' @importFrom stats weighted.mean
#' @importFrom scales percent
#' @rdname plotDupLevels-methods
#' @export
setMethod(
    "plotDupLevels",signature = "FastpDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        pwfCols, warn = 20, fail = 50, showPwf = FALSE, plotlyLegend = FALSE,
        plotType = c("bar", "heatmap"), barFill = "blue", barCol = "blue",
        cluster = FALSE, dendrogram = FALSE,  scaleFill = NULL,
        plotTheme = theme_get(), heat_w = 8, maxLevel = 10, ...
    ){

        mod <- "Duplication"
        df <- getModule(x, mod)
        plotType <- match.arg(plotType)

        if (!length(df)) {
            p <- .emptyPlot("No Duplication Levels Detected")
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(df, labels, pattern = pattern)
        labels <- labels[names(labels) %in% df$Filename]
        key <- names(labels)

        ## Sort out the colours
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))

        ## Given that levels are calculated across R1/2 and are simply summarised
        ## build the plot as a simply barplot
        if (plotType == "bar") {

            df$Filename <- forcats::fct_inorder(labels[df$Filename])
            df$x <- as.integer(df$Filename)
            df$Rate <- percent(df$rate, 0.01)
            pwfCols <- setAlpha(pwfCols, 0.2)

            ## Set the limits & rectangles
            ylim <- c(0, max(df$rate) * 1.05)
            expand_x <- round(0.015*(max(df$x) - min(df$x)), 1)
            rects <- tibble(
                xmin = min(df$x - 0.5) - expand_x,
                xmax = max(df$x + 0.5) + expand_x,
                ymin = c(0, warn, fail) / 100,
                ymax = c(warn, fail, max(ylim * 100)) / 100,
                Status = forcats::fct_inorder(c("PASS", "WARN", "FAIL"))
            )
            xlim <- c(rects$xmin[[1]], rects$xmax[[1]])
            rects <- dplyr::filter(rects, ymin < max(ylim))
            rects$ymax[rects$ymax > max(ylim)] <- max(ylim)

            p <- ggplot(df, aes(label = Filename, rate = Rate))
            if (showPwf) {
                p <- p + geom_rect(
                    data = rects,
                    aes(
                        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        fill = Status
                    ),
                    inherit.aes = FALSE
                )
            }

            p <- p +
                geom_col(
                    aes(x, !!sym("rate")), fill = barFill, colour = barCol
                ) +
                scale_x_continuous(
                    breaks = seq_along(labels), labels = df$Filename,
                    expand = c(0, 0),limits = xlim
                ) +
                scale_y_continuous(
                    expand = expansion(c(0, 0)), labels = percent, limits = ylim
                ) +
                scale_fill_manual(values = getColours(pwfCols)) +
                theme_bw() +
                ggtitle("Duplication Rate") +
                labs(x = "Filename", y = "Duplication Rate (%)")
            p <- .updateThemeFromDots(p, ...)
            if (usePlotly) {
                hv <- c("Filename", "Rate")
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                p <- plotly::ggplotly(p, tooltip = hv)
                p$x$data <- lapply(p$x$data, .hidePWFRects)
            }
        }

        if (plotType == "heatmap") {

            ## Extract the histogram data
            df <- unnest(df, !!sym("histogram"))
            if (!nrow(df)) {
                msg <- "No Duplication Histogram Data Found"
                message(msg)
                p <- .emptyPlot(msg)
                if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
                return(p)
            }

            ## Prepare for plotting
            df <- dplyr::arrange(df, !!sym("duplication_level"))
            df$duplication_level <- ifelse(
                df$duplication_level > maxLevel,
                paste0(">", maxLevel), as.character(df$duplication_level)
            )
            df$duplication_level <- forcats::fct_inorder(df$duplication_level)
            df <- group_by(df, Filename, !!sym("duplication_level"))
            df <- summarise(
                df,
                rate = unique(!!sym("rate")),
                histogram = sum(!!sym("histogram")),
                mean_gc = weighted.mean(
                    !!sym("mean_gc"), w = !!sym("duplication_rate")
                ),
                duplication_rate = sum(!!sym("duplication_rate")),
                .groups = "drop_last"
            )
            df <- mutate(
                df,
                xmax = cumsum(!!sym("duplication_rate")),
                xmin = c(0, xmax[-dplyr::n()])
            )
            df <- ungroup(df)

            ## Set up the dendrogram
            n <- length(unique(df$Filename))
            dx <- list()
            dx$segments <- tibble(
                x = numeric(), y = numeric(), xend = numeric(), yend = numeric()
            )
            if (n > 1) {
                ## This will only run if we have > 1 sample
                clusterDend <- .makeDendro(
                    df, "Filename", "duplication_level", "duplication_rate"
                )
                dx <- ggdendro::dendro_data(clusterDend)
                if (dendrogram | cluster) key <- labels(clusterDend)
                if (!dendrogram) dx$segments <- dx$segments[0,]
            } else {
                dendrogram <- FALSE
            }
            df$Filename <- factor(labels[df$Filename], levels = labels[key])

            ## Now add the correct y-coordinates
            df$y <- as.integer(df$Filename) - 0.5
            df$ymin <- df$y - 0.5
            df$ymax <- df$y + 0.5

            ## The PWF status
            status_df <- dplyr::distinct(df, Filename, !!sym("rate"))
            status_df$Status <- cut(
                status_df$rate, breaks = c(0, warn, fail, 100) / 100,
                labels = c("PASS", "WARN", "FAIL")
            )
            if (!showPwf) status_df <- status_df[0, ]

            if (is.null(scaleFill))
                scaleFill <- scale_fill_viridis_d(option = "inferno")
            stopifnot(is(scaleFill, "ScaleDiscrete"))
            stopifnot(scaleFill$aesthetics == "fill")
            stopifnot(is(plotTheme, "theme"))

            ## Tidy up for plotting
            dl <- "Duplication Level"
            names(df)[names(df) == "duplication_level"] <- dl
            dupLevels <- levels(df[[dl]])
            df$Rate <- percent(df$duplication_rate, 0.01)
            hj <-  0.5 * heat_w / (heat_w + 1 + dendrogram)
            p <- ggplot(df, aes(rate = Rate, name = Filename)) +
                geom_rect(
                    aes(
                        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        fill = !!sym(dl)
                    )
                ) +
                scale_x_continuous(expand = expansion(0, 0), labels = percent) +
                scale_y_continuous(
                    breaks = unique(df$y), labels = levels(df$Filename),
                    expand = c(0, 0),position = "right"
                ) +
                scaleFill +
                labs(
                    x = "% of Library", y = "Filename",
                    fill = "Duplication\nLevel"
                ) +
                plotTheme

            hv <- c("Filename", dl, "Rate")
            ## Calling plotly with 1 sample goes a bit funny. Not sure why
            if (usePlotly) p <- p + labs(y = NULL)
            p <- .prepHeatmap(
                p, status_df, dx$segments, usePlotly, heat_w, pwfCols, hv
            )

        }
        p
    }
)
