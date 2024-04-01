#' @title Draw an N Content Plot
#'
#' @description Draw an N Content Plot across one or more FastQC reports
#'
#' @details
#' This extracts the N_Content from the supplied object and generates a ggplot2
#' object, with a set of minimal defaults.
#' The output of this function can be further modified using the standard
#' ggplot2 methods.
#'
#' When `x` is a single FastqcData object line plots will always be drawn
#' for all Ns.
#' Otherwise, users can select line plots or heatmaps.
#'
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical`. Output as ggplot2 (default) or plotly
#' object.
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' @param pattern Regex used to trim the end of filenames
#' @param warn,fail The default values for warn and fail are 5 and 10
#' respectively (i.e. percentages)
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param showPwf logical(1) Show the PASS/WARN/FAIL status
#' @param module Used for Fastp* structures to show results before or after
#' filtering
#' @param reads Show plots for read1, read2 or both.
#' @param moduleBy,readsBy How to show each module or set of reads on the plot
#' @param lineCol Defaults to red
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param heat_w Relative width of any heatmap plot components
#' @param scaleFill,scaleColour,scaleLine ggplot2 scale objects
#' @param plotTheme \link[ggplot2]{theme} object
#' @param plotlyLegend logical(1) Show legend on interactive plots
#' @param lineCol Line colours
#' @param ... Used to pass additional attributes to theme() for FastqcData
#' objects and to geom* calls for FastpData-based objects
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @examples
#'
#' ## Using a Fastp Data object
#' fl <- system.file("extdata/fastp.json.gz", package = "ngsReports")
#' fp <- FastpData(fl)
#' plotNContent(fp)
#' plotNContent(
#'   fp, pattern = "_001.+",
#'   moduleBy = "colour", scaleColour = scale_colour_brewer(palette = "Set1"),
#'   plotTheme = theme(
#'     legend.position = 'inside', legend.position.inside = c(0.99, 0.99),
#'     legend.justification = c(1, 1), plot.title = element_text(hjust = 0.5)
#'   )
#' )
#'
#' @docType methods
#'
#' @import ggplot2
#' @importFrom dplyr vars
#' @import tibble
#' @importFrom tidyr unnest
#' @importFrom stringr str_split
#'
#' @name plotNContent
#' @rdname plotNContent-methods
#' @export
setGeneric(
    "plotNContent",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...) standardGeneric("plotNContent")
)
#' @rdname plotNContent-methods
#' @export
setMethod(
    "plotNContent", signature = "ANY",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){.errNotImp(x)}
)
#' @rdname plotNContent-methods
#' @export
setMethod(
    "plotNContent", signature = "FastqcData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
        warn = 5, fail = 20, showPwf = TRUE, ..., lineCol = "red"
    ){

        ## Get the NContent
        df <- getModule(x, "Per_base_N_content")
        colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

        ## Handle empty/missing modules
        msg <- c()
        if (!length(df)) msg <- "No N Content Module Detected"
        if (sum(df[["Percentage"]]) == 0) msg <- "No N Content in Sequences"
        if (!is.null(msg)) {
            p <- ngsReports:::.emptyPlot(msg)
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        ## Sort out the colours
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))
        pwfCols <- setAlpha(pwfCols, 0.2)

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern)
        labels <- labels[names(labels) %in% df$Filename]
        df$Filename <- labels[df$Filename]
        df$Base <- factor(df$Base, levels = unique(df$Base))
        df$xValue <- as.integer(df$Base)

        ## Setup the BG colours
        rects <- tibble(
            xmin = 0, xmax = max(df$xValue),
            ymin = c(0, warn, fail), ymax = c(warn, fail, 100),
            Status = c("PASS", "WARN", "FAIL")
        )

        yLab <- "N Content (%)"
        x <- "xValue"
        p <- ggplot(df)
        if (showPwf)
            p <- p + geom_rect(
                data = rects,
                aes(
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = Status
                )
            )

        p <- p +
            geom_line(aes(!!sym(x), Percentage), colour = lineCol) +
            geom_point(
                aes(!!sym(x), Percentage, group = Base),
                size = 0, colour = rgb(0, 0, 0, 0)
            ) +
            scale_fill_manual(values = getColours(pwfCols)) +
            scale_x_continuous(
                breaks = unique(df$xValue), labels = levels(df$Base),
                expand = c(0,0)
            ) +
            scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
            facet_wrap(~Filename) +
            labs(x = "Position in Read (bp)", y = yLab) +
            guides(fill = "none") +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )

        ## Add the basic customisations
        p <- .updateThemeFromDots(p, ...)

        if (usePlotly) {
            p <- p +
                xlab("") +
                theme(legend.position = "none")
            p <- suppressMessages(plotly::ggplotly(p))
            p <- plotly::layout(p, yaxis1 = list(title = yLab))

            ## Set the hoverinfo for bg rectangles to the vertices only,
            ## This will effectively hide them
            p$x$data <- lapply(p$x$data, .hidePWFRects)
            ## Hide the xValue parameter to make it look nicer
            p$x$data[[6]]$text <- gsub(
                "(.+)(xValue.+)(Percentage.+)", "\\1\\3", p$x$data[[6]]$text
            )
        }
        p
    }
)
#' @rdname plotNContent-methods
#' @export
setMethod(
    "plotNContent", signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
        warn = 5, fail = 20, showPwf = TRUE, cluster = FALSE,
        dendrogram = FALSE, heat_w = 8, scaleFill = NULL, ...
    ){

        ## Get the NContent
        df <- getModule(x, "Per_base_N_content")
        colnames(df) <- gsub("N-Count", "Percentage", colnames(df))

        ## Handle empty/missing modules
        msg <- c()
        if (!length(df)) msg <- "No N Content Module Detected"
        if (sum(df[["Percentage"]]) == 0) msg <- "No N Content in Sequences"
        if (!is.null(msg)) {
            p <- .emptyPlot(msg)
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        ## Sort out the colours
        if (missing(pwfCols)) pwfCols <- ngsReports::pwf
        stopifnot(.isValidPwf(pwfCols))

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern)
        labels <- labels[names(labels) %in% df$Filename]

        ## fill bins up to the max sequence length
        df$Base <- lapply(
            df$Base,
            function(x){
                rng <- as.integer(str_split(x, pattern = "-")[[1]])
                seq(min(rng), max(rng), by = 1L)
            }
        )
        df <- unnest(df, Base)

        ## Now define the order for a dendrogram if required
        key <- names(labels)
        cols <- c("Filename", "Base", "Percentage")
        clusterDend <- .makeDendro(df[cols], "Filename", "Base", "Percentage")
        dx <- ggdendro::dendro_data(clusterDend)
        if (dendrogram | cluster) key <- labels(clusterDend)
        if (!dendrogram) dx$segments <- dx$segments[0,]
        df$Filename <- factor(labels[df$Filename], levels = labels[key])

        if (is.null(scaleFill)) {
            cols <- .makePwfGradient(
                df$Percentage, pwfCols,
                breaks = c(0, warn, fail, 101), passLow = TRUE,
                na.value = "white"
            )
            scaleFill <- do.call("scale_fill_gradientn", cols)
        }
        stopifnot(is(scaleFill, "ScaleContinuous"))
        stopifnot(scaleFill$aesthetics == "fill")

        xLab <- "Position in Read (bp)"
        hj <- 0.5 * heat_w / (heat_w + 1 + dendrogram)
        p <- ggplot(df, aes(Base, Filename, fill = Percentage, label = Base)) +
            geom_tile() +
            scaleFill +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0), position = "right") +
            labs(x = xLab, y = NULL, fill = "%N") +
            ggtitle("Per Base N Content") +
            theme_bw() +
            theme(
                panel.background = element_blank(),
                plot.title = element_text(hjust = hj),
                axis.title.x = element_text(hjust = hj)
            )
        if (dendrogram | showPwf)
            p <- p + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))
        p <- .updateThemeFromDots(p, ...)

        ## Reset the status using current values
        status <- dplyr::summarise(
            dplyr::group_by(df, Filename),
            Percentage = max(Percentage, na.rm = TRUE)
        )
        status$Status <- cut(
            status$Percentage, include.lowest = TRUE,
            breaks = c(0, warn, fail, 101), labels = c("PASS", "WARN", "FAIL")
        )
        status <- subset(status, Filename %in% key)
        status$Filename <- factor(labels[status$Filename], levels = labels[key])
        if (!showPwf) status <- status[0,]

        .prepHeatmap(p, status, dx$segments, usePlotly, heat_w, pwfCols)

    }
)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr unnest
#' @importFrom tidyselect any_of
#' @importFrom scales label_percent comma
#' @importFrom rlang sym !!
#' @rdname plotNContent-methods
#' @export
setMethod(
    "plotNContent", signature = "FastpData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        module = c("Before_filtering", "After_filtering"),
        moduleBy = c("facet", "colour", "linetype"),
        reads = c("read1", "read2"), readsBy = c("facet", "colour", "linetype"),
        scaleColour = NULL, scaleLine = NULL, plotTheme = theme_get(),
        plotlyLegend = FALSE, ...
    ){

        ## Check args and load
        mod <- match.arg(module, several.ok = TRUE) # We can plot B4/After
        reads <- match.arg(reads, several.ok = TRUE)
        data <- lapply(mod, function(m) getModule(x, m))
        names(data) <- mod
        data <- lapply(data, bind_rows, .id = "reads")
        df <- bind_rows(data,.id = "Module")
        cols <- c(
            "Module", "reads", "Filename", "fqName", "total_reads",
            "content_curves", "N", "position"
        )
        df <- dplyr::select(df, any_of(cols))
        df <- unnest(df, !!sym("content_curves"))
        df <- dplyr::select(df, any_of(cols))

        ## Make a blank plot if no data is found
        if (!length(df)) {
            msg <- "No N Content in Reports"
            p <- .emptyPlot(msg)
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        ## Drop the suffix, or check the alternate labels
        lbl_df <- dplyr::distinct(df, Filename, fqName)
        fnames <- .makeLabels(lbl_df, labels, pattern = pattern)
        df$Filename <- factor(fnames[df$Filename], levels = fnames)
        df$Filename <- fct_rev(df$Filename)
        fqNames <- .makeLabels(
            lbl_df, labels, pattern = pattern, col = "fqName"
        )
        df$fqName <- fqNames[df$fqName]
        df$Module <- factor(df$Module, levels = mod)
        df$module <- df$Module ## Makes it easier for users to make typos

        ## Sort out plotting args
        readsBy <- match.arg(readsBy)
        moduleBy <- match.arg(moduleBy)
        if (readsBy == moduleBy & readsBy != "facet") stop(
            "Cannot set the same plotting parameter to both reads and module"
        )
        lt <- col <- NULL
        if (readsBy == "linetype") lt <- sym("fqName")
        if (moduleBy == "linetype") lt <- sym("Module")
        if (readsBy == "colour") col <- sym("fqName")
        if (moduleBy == "colour") col <- sym("Module")
        fm <- dplyr::case_when(
            readsBy == "facet" & moduleBy == "facet" ~ "Module ~ fqName",
            readsBy != "facet" & moduleBy == "facet" ~ ". ~ Module",
            readsBy == "facet" & moduleBy != "facet" ~ ". ~ fqName",
            TRUE ~ ". ~ ."
        )
        fm <- as.formula(fm)
        if (!is.null(scaleColour)) {
            stopifnot(is(scaleColour, "ScaleDiscrete"))
            stopifnot(scaleColour$aesthetics == "colour")
        }
        if (!is.null(scaleLine)) {
            stopifnot(is(scaleLine, "ScaleDiscrete"))
            stopifnot(scaleLine$aesthetics == "linetype")
        }
        stopifnot(is(plotTheme, "theme"))

        ## Final mods for plotting
        df[["% N"]] <- round(100 * df$N, 2)
        df[["N Reads"]] <- comma(df$N * df$total_reads, 1)
        names(df) <- gsub("position", "Position", names(df))
        df_cols <- colnames(df)
        p <- ggplot(df, aes(label = !!sym("N Reads"))) +
            geom_line(
                aes(
                    Position, !!sym("% N"), colour = {{ col }},
                    linetype = {{ lt }}
                ), ...
            ) +
            ggtitle("N Content") +
            scaleColour + scaleLine +
            scale_y_continuous(
                labels = label_percent(scale = 1),
                expand = expansion(c(0.01, 0.05))
            ) +
            facet_grid(fm) +
            labs(y = "% N") +
            theme_bw() +
            plotTheme

        if (usePlotly) {
            if (!plotlyLegend) p <- p + theme(legend.position = "none")
            hv <- c("fqName", "Position", "% N", "label", "colour", "linetype")
            p <- suppressWarnings(plotly::ggplotly(p, tooltip = hv))
        }
        p
    }
)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr unnest
#' @importFrom tidyselect any_of
#' @importFrom scales percent comma
#' @importFrom rlang sym !!
#' @rdname plotNContent-methods
#' @export
setMethod(
    "plotNContent", signature = "FastpDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        module = c("Before_filtering", "After_filtering"),
        reads = c("read1", "read2"), scaleFill = NULL, plotTheme = theme_get(),
        cluster = FALSE, dendrogram = FALSE, heat_w = 8, ...
    ){

        ## Check args
        mod <- match.arg(module, several.ok = TRUE)
        # We can't plot B4/After if we cluster
        reads <- match.arg(reads, several.ok = TRUE)
        n_facets <- length(mod) * length(reads)
        if ((dendrogram | cluster) & n_facets > 2) {
            message(
                "Cannot cluster when plotting both modules and both sets of reads. ",
                "Setting cluster and dendrogram to FALSE"
            )
            cluster <- dendrogram <- FALSE
        }

        ## Setup the data
        data <- lapply(mod, function(m) getModule(x, m)[reads])
        names(data) <- mod
        data <- lapply(data, bind_rows, .id = "reads")
        df <- bind_rows(data, .id = "Module")
        cols <- c(
            "Module", "reads", "Filename", "fqName", "total_reads",
            "content_curves", "N", "position"
        )
        df <- dplyr::select(df, any_of(cols))
        df <- unnest(df, !!sym("content_curves"))
        df <- dplyr::select(df, any_of(cols))

        ## Make a blank plot if no data is found
        if (!length(df)) {
            msg <- "No N Content in Reports"
            p <- .emptyPlot(msg)
            if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
            return(p)
        }

        ## Drop the suffix, or check the alternate labels
        lbl_df <- dplyr::distinct(df, Filename, fqName)
        fqNames <- .makeLabels(
            lbl_df, labels, pattern = pattern, col = "fqName"
        )
        df$fqName <- fqNames[df$fqName]
        labels <- .makeLabels(lbl_df, labels, pattern = pattern)
        key <- names(labels)
        labels <- labels[key %in% df$Filename]

        ## Set up the dendrogram & labels
        n <- length(labels)
        if (n == 1 & (cluster | dendrogram))
            message("Cannot cluster one file. Ignoring cluster and dendgrogram")
        if (n > 1 & (cluster | dendrogram)) {
            df$rowVal <- paste(df$position, df$reads, df$Module)
            clusterDend <- .makeDendro(df, "Filename", "rowVal", "N")
            dx <- ggdendro::dendro_data(clusterDend)
            if (dendrogram | cluster) key <- labels(clusterDend)
        } else {
            cluster <- dendrogram <- FALSE
            dx <- list()
            dx$segments <- lapply(rep_len(0, 4), numeric)
            names(dx$segments) <- c("x", "y", "xend", "yend")
            dx$segments <- as_tibble(dx$segments)
        }
        if (!dendrogram) dx$segments <- dx$segments[0,]
        ## Now set everything as factors
        lv <- labels[key]
        df$Filename <- factor(labels[df$Filename], levels = lv)

        if (is.null(scaleFill)) scaleFill <- scale_fill_viridis_c()
        stopifnot(is(scaleFill, "ScaleContinuous"))
        stopifnot(scaleFill$aesthetics == "fill")
        stopifnot(is(plotTheme, "theme"))

        ## Tidy up for plotting
        df$Module <- factor(df$Module, levels = mod)
        df[["% N"]] <- percent(df$N, 0.01)
        df[["N Reads"]] <- comma(df$N * df$total_reads, 1)
        names(df) <- gsub("position", "Position", names(df))
        fm <- as.formula(". ~ reads")
        if (n_facets == 2 & length(mod) == 2) fm <- as.formula(". ~ Module")
        if (n_facets == 4) fm <- as.formula("Module ~ reads")
        p <- ggplot(
            df,
            aes(
                Position, Filename, fill = !!sym("N"), label = !!sym("fqName"),
                percent = !!sym("% N"), total = !!sym("N Reads")
            )
        ) +
            geom_raster(...) +
            facet_grid(fm, switch = "y") +
            ggtitle("N Content") +
            scaleFill +
            scale_y_discrete(expand = rep(0, 4), position = "right") +
            scale_x_continuous(expand = rep(0, 4)) +
            theme_bw() +
            plotTheme
        if (dendrogram)
            p <- p + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))

        if (usePlotly) p <- p + theme(
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
        tt <- c("fqName", "Position", "% N", "N Reads", "Module")
        p <- .prepHeatmap(
            p, tibble(), usePlotly, heat_w, segments = dx$segments, hv = tt
        )
        p
    }
)
