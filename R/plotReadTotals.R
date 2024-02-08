#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' The raw data from [readTotals()] can otherwise be used to manually
#' create a plot.
#'
#' Duplication levels are based on the value shown on FASTQC reports at the
#' top of the DeDuplicatedTotals plot, which is known to be inaccurate.
#' As it still gives a good guide as to sequence diversity it is included as
#' the default. This can be turned off by setting `duplicated = FALSE`.
#'
#' For FastpDataList objects, duplication statistics are not part of the
#' default module containing ReadTotals. However, the status of reads and the
#' reason for being retained or filtered is, and as such these are shown
#' instead of duplication statistics.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical` Default `FALSE` will render using
#' ggplot. If `TRUE` plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' @param pattern Regex used to trim the end of filenames
#' @param duplicated logical(1). Include deduplicated read total estimates
#' to plot charts
#' @param status logical(1) Include read status in the plot
#' @param bars If `duplicated = TRUE`, show unique and deduplicated reads
#' as "stacked" or "adjacent".
#' @param barCols Colours for duplicated and unique reads.
#' @param expand.y Passed to [`ggplot2::expansion`] for the axis showing totals
#' @param vertBars logical(1) Show bars as vertical or horizontal
#' @param adjPaired Scale read totals by 0.5 when paired
#' @param divBy Scale read totals by this value. The default shows the y-axis
#' in millions for FastpDataList objects, but does not scale FastQC objects,
#' for the sake of backwards compatability
#' @param scaleFill ScaleDiscrete function to be applied to the plot
#' @param labMin Only show labels for filtering categories higher than this
#' values as a proportion of reads. Set to any number > 1 to turn off labels
#' @param labelVJ Relative vertical position to labels within each bar.
#' @param labelFill Passed to \link[ggplot2]{geom_label}
#' @param plotTheme \link[ggplot2]{theme} to be added to the plot
#' @param plotlyLegend logical(1) Show legend on interactive plots
#' @param ... Used to pass additional attributes to theme()
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
#' # Plot the Read Totals showing estimated duplicates
#' plotReadTotals(fdl)
#'
#' # Plot the Read Totals without estimated duplicates
#' plotReadTotals(fdl, duplicated = FALSE)
#'
#' @return Returns a ggplot or plotly object
#'
#' @docType methods
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom tidyselect one_of
#' @importFrom forcats fct_inorder
#'
#' @name plotReadTotals
#' @rdname plotReadTotals-methods
#' @export
setGeneric(
    "plotReadTotals",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){
        standardGeneric("plotReadTotals")
    }
)
#' @rdname plotReadTotals-methods
#' @export
setMethod(
    "plotReadTotals", signature = "ANY",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){
        .errNotImp(x)
    })
#' @importFrom scales comma
#' @rdname plotReadTotals-methods
#' @export
setMethod(
    "plotReadTotals", signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels,  pattern = ".(fast|fq|bam).*",
        duplicated = TRUE, bars = c("stacked", "adjacent"), vertBars = TRUE,
        divBy = 1, barCols = c("red","blue"), expand.y = c(0, 0.02),
        plotlyLegend = FALSE, ...
    ){

        stopifnot(is.logical(duplicated))
        df <- readTotals(x)

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern = pattern, ...)
        df$Filename <- labels[df$Filename]
        df$Filename <- forcats::fct_inorder(df$Filename)

        ## Get the colours for the barplot
        barCols <- rep_len(barCols, duplicated + 1)
        y_col <- "Read Totals"
        if (divBy != 1) y_col <- paste0("Read Totals (x", comma(divBy, 1), ")")
        df[[y_col]] <- df$Total_Sequences / divBy

        ## Check the axis expansion
        stopifnot(is.numeric(expand.y))
        expand.y <- expansion(rep_len(expand.y, 2))
        yMax <- max(df[[y_col]])

        if (!duplicated) {

            p <- ggplot(df, aes(x = Filename, y = !!sym(y_col))) +
                geom_bar(stat = "identity", fill = barCols) +
                scale_y_continuous(
                    labels = comma, limits = c(0, yMax), expand = expand.y
                ) +
                theme_bw()
            if (vertBars) p <- p + coord_flip()
            p <- .updateThemeFromDots(p, ...)
            if (usePlotly) p <- plotly::ggplotly(p)

        }

        if (duplicated) {

            ## Add the information to a joined data.frame
            deDup <- getModule(x, "Total_Deduplicated_Percentage")
            deDup$Filename <- labels[deDup$Filename]
            deDup$Filename <- factor(
                deDup$Filename, levels = levels(df$Filename)
            )
            names(deDup) <- gsub("Total_Deduplicated_", "", names(deDup))
            df <- dplyr::left_join(deDup, df, by = "Filename")

            ##Setup the df for plotting
            types <- c("Unique", "Duplicated")
            df$Unique <- df$Percentage * df[[y_col]] / 100
            df$Unique <- round(df$Unique, log10(divBy))
            df$Duplicated <- df[[y_col]] - df$Unique
            df <- df[c("Filename", types)]
            df <- tidyr::gather(df, "Type", !!sym(y_col), one_of(types))

            bars <- match.arg(bars)
            barPos <- c(adjacent = "dodge", stacked = "stack")[bars]

            ## Make the plot
            p <- ggplot(df, aes(x = Filename, y = !!sym(y_col), fill = Type)) +
                geom_bar(stat = "identity", position = barPos) +
                scale_y_continuous(
                    labels = comma, limits = c(0, yMax), expand = expand.y
                ) +
                scale_fill_manual(values = barCols) +
                theme_bw()
            if (vertBars) p <- p + coord_flip()
            p <- .updateThemeFromDots(p, ...)

            if (usePlotly) {

                # Hide the legend
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                # Render as a plotly object
                p <- suppressMessages(suppressWarnings(plotly::ggplotly(p)))
            }

        }

        ## Draw the plot
        p
    }
)
#' @importFrom dplyr mutate
#' @importFrom scales percent comma
#' @rdname plotReadTotals-methods
#' @export
setMethod(
    "plotReadTotals", signature = "FastpDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        adjPaired  = TRUE, divBy = 1e6, scaleFill = NULL, labMin = 0.05,
        status = TRUE, labelVJ = 0.5, labelFill = "white",
        plotTheme = theme_get(), vertBars = FALSE, plotlyLegend = FALSE,
        expand.y = c(0, 0.05), ...
    ){

        ## Sort out the summary data
        module <- "Summary"
        data <- getModule(x, module)
        df <- data$Filtering_result
        df <- dplyr::filter(df, !!sym("total") > 0)
        if (!nrow(df)) {
            ## If no filtering was performed, just use the Before_filtering data
            df <- data$Before_filtering
            stopifnot(nrow(df))
            df$total <- df$total_reads
            df$result <- "total_reads"
            df$rate <- 1
        }

        ## Get the correct order for filling columns
        levels <- sort(rank(df$total))
        levels <- unique(names(levels))
        df$result <- factor(df$result, levels = levels)
        df$result <- forcats::fct_relabel(
            df$result,
            \(x) str_remove_all(str_to_title(gsub("_", " ", x)), " Reads")
        )
        paired <- getModule(x, "paired")
        df <- dplyr::left_join(df, paired, by = "Filename")
        y_col <- "Read Totals"
        if (divBy != 1) y_col <- paste0("Read Totals (x", comma(divBy, 1), ")")
        df[[y_col]] <- df$total / divBy
        if (adjPaired) df[[y_col]][df$paired] <- 0.5 * df[[y_col]][df$paired]
        df <- mutate(df, cumsum = cumsum(!!sym(y_col)), .by = Filename)
        df$label_y <- df$cumsum - (1 - labelVJ) * df[[y_col]]

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(df, labels, pattern = pattern, ...)
        df$Filename <- factor(labels[df$Filename], levels = labels)

        ## Add additional columns for a nice plotly figure
        df$Total <- df$total
        if (adjPaired) df$Total[df$paired] <- 0.5 * df$Total[df$paired]
        df$Total <- comma(df$Total)
        df[["% Reads"]] <- percent(df$rate, 0.01)
        fill_col <- "Filtering Result"
        names(df) <- gsub("result", fill_col, names(df))
        if (is.null(scaleFill)) scaleFill <- scale_fill_viridis_d(
            option = "cividis", direction = -1
        )
        stopifnot(is(scaleFill, "ScaleDiscrete"))
        stopifnot(scaleFill$aesthetics == "fill")
        stopifnot(is(plotTheme, "theme"))
        fill_col <- sym(fill_col)
        if (!status) fill_col <- NULL

        p <- ggplot(
            df,
            aes(
                Filename, !!sym(y_col), fill = {{ fill_col }},
                total = !!sym("Total"),
                paired = !!sym("paired"), percent = !!sym("% Reads")
            )
        ) +
            geom_col(...) +
            scale_y_continuous(
                labels = comma, expand = expansion(rep_len(expand.y, 2))
            ) +
            scaleFill +
            theme_bw() +
            plotTheme
        if (vertBars) p <- p + coord_flip()

        if (!usePlotly) {
            if (status) p <- p + geom_label(
                aes(
                    Filename, y = !!sym("label_y"),
                    label = percent(!!sym("rate"), 0.1)
                ),
                data = dplyr::filter(df, !!sym("rate") >= labMin),
                fill = labelFill,
                ...
            )
        } else {
            if (!plotlyLegend) p <- p + theme(legend.position = "none")
            tt <- c("Filename", "Filtering Result", "Total", "% Reads")
            p <- plotly::ggplotly(p, tooltip = tt)
        }

        p
    }
)
