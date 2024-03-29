#' @title Plot the per base content as a heatmap
#'
#' @description Plot the Per Base content for a set of FASTQC files.
#'
#' @details
#' Per base sequence content (%A, %T, %G, %C), is shown as four overlaid
#' heatmap colours when plotting from multiple reports. The individual line
#' plots are able to be generated by setting `plotType = "line"`, and the
#' layout is determined by `facet_wrap` from ggplot2.
#'
#' Individual line plots are also generated when plotting from a single
#' `FastqcData` object.
#'
#' If setting `usePlotly = TRUE` for a large number of reports, the plot
#' can be slow to render.
#' An alternative may be to produce a plot of residuals for each base, produced
#' by taking the position-specific mean for each base.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#' @param usePlotly `logical`. Generate an interactive plot using plotly
#' @param labels An optional named vector of labels for the file names.
#' All file names must be present in the names of the vector.
#' @param pattern Regex to remove from the end of any filenames
#' @param plotType `character`. Type of plot to generate. Must be "line",
#' "heatmap" or "residuals"
#' @param pwfCols Object of class [PwfCols()] to give colours for
#' pass, warning, and fail
#' values in plot
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param heat_w Relative width of any heatmap plot components
#' @param module Fastp Module to show. Can only be Before/After_filtering
#' @param reads Which set of reads to show
#' @param readsBy,moduleBy When plotting both R1 & R2 and both modules,
#' separate by either linetype or linetype
#' @param bases Which bases to draw on the plot. Also becomes the default
#' plotting order by setting these as factor levels
#' @param scaleColour Discrete colour scale as a ggplot ScaleDiscrete object
#' If not provided, will default to \link[ggplot2]{scale_colour_manual}
#' @param scaleLine Discrete scale_linetype object. Only relevant if plotting
#' values by linetype
#' @param plotTheme \link[ggplot2]{theme} object to be applied. Note that all
#' plots will have \link[ggplot2]{theme_bw} theme applied by default, as well as
#' any additional themes supplied here
#' @param plotlyLegend logical(1) Show legends for interactive plots. Ignored
#' for heatmaps
#' @param showPwf Show PASS/WARN/FAIL categories as would be defined in a FastQC
#' report
#' @param expand.x,expand.y Passed to \link[ggplot2]{expansion} in the x- and
#' y-axis scales respectively
#' @param ... Used to pass additional attributes to plotting geoms
#' @param nc Specify the number of columns if plotting a FastqcDataList as line
#' plots. Passed to \link[ggplot2]{facet_wrap}.
#' @param warn,fail Default values for WARN and FAIL based on FastQC reports.
#' Only applied to heatmaps for FastpDataList objects
#'
#' @return A ggplot2 object or an interactive plotly object
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
#' plotSeqContent(fdl)
#'
#' fp <- FastpData(system.file("extdata/fastp.json.gz", package = "ngsReports"))
#' plotSeqContent(fp)
#' plotSeqContent(fp, moduleBy = "linetype", bases = c("A", "C", "G", "T"))
#'
#' @docType methods
#'
#' @importFrom grDevices rgb
#' @importFrom dplyr mutate vars group_by ungroup left_join
#' @importFrom scales percent
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect one_of all_of
#' @import ggplot2
#'
#' @name plotSeqContent
#' @rdname plotSeqContent-methods
#' @export
setGeneric(
    "plotSeqContent",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...) standardGeneric("plotSeqContent")
)
#' @rdname plotSeqContent-methods
#' @export
setMethod(
    "plotSeqContent", signature = "ANY",
    function(x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", ...){.errNotImp(x)}
)
#' @importFrom scales label_percent
#' @rdname plotSeqContent-methods
#' @export
setMethod(
    "plotSeqContent", signature = "FastqcData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        bases = c("A", "T", "C", "G"), scaleColour = NULL,
        plotTheme = theme_get(), plotlyLegend = FALSE, expand.x = 0.02,
        expand.y = c(0, 0.05), ...
    ){

        ## Get the SequenceContent
        df <- getModule(x, "Per_base_sequence_content")
        names(df)[names(df) == "Base"] <- "Position"

        if (!length(df)) {
            p <- .emptyPlot("No Sequence Content Module Detected")
            if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
            return(p)
        }

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern = pattern, ...)
        labels <- labels[names(labels) %in% df$Filename]
        df$Filename <- labels[df$Filename]

        bases <- match.arg(bases, several.ok = TRUE)
        df <- tidyr::pivot_longer(
            df, cols = all_of(bases), names_to = "Base", values_to = "Percent"
        )
        df$Base <- factor(df$Base, levels = bases)
        df$Position <- as.integer(str_extract(df$Position, "^[0-9]+"))
        df$Percent <- round(df$Percent, 2)
        ## This was a good idea but looks ugly. Maybe make it an option??
        # df <- tidyr::complete(
        #   df, Position = seq_len(max(df$Position)), nesting(Filename, Base)
        # )
        # df <- dplyr::arrange(df, Filename, Base, Position)
        # df <- zoo::na.locf(df)

        ## set colours & theme
        if (is.null(scaleColour)) {
            baseCols <- c(
                `T` = "red", G = "black", A = "green", C = "blue"
            )[bases]
            scaleColour <- scale_colour_manual(values = baseCols)
        }
        stopifnot(is(scaleColour, "ScaleDiscrete"))
        stopifnot(scaleColour$aesthetics == "colour")
        stopifnot(is(plotTheme, "theme"))

        xLab <- "Position in read (bp)"
        yLab <- "Percent"
        p <- ggplot(
            df, aes(Position, Percent, label = Position, colour = Base)
        ) +
            geom_line(...) +
            facet_wrap(~Filename) +
            scale_y_continuous(
                limits = c(0, max(df$Percent)),
                labels = label_percent(scale = 1),
                expand = expansion(rep_len(expand.y, 2))
            ) +
            scale_x_continuous(expand = expansion(rep_len(expand.x, 2))) +
            scaleColour +
            guides(fill = "none") +
            labs(x = xLab, y = yLab) +
            theme_bw() + plotTheme

        if (usePlotly) {
            ttip <- c("y", "label", "colour")
            if (!plotlyLegend) p <- p + theme(legend.position = "none")
            p <- plotly::ggplotly(p, tooltip = ttip)
        }

        p

    }
)
#' @rdname plotSeqContent-methods
#' @export
setMethod(
    "plotSeqContent", signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*", pwfCols,
        showPwf = TRUE, plotType = c("heatmap", "line", "residuals"),
        scaleColour = NULL, plotTheme = theme_get(), cluster = FALSE,
        dendrogram = FALSE, heat_w = 8, plotlyLegend = FALSE, nc = 2, ...
    ){

        ## Get the SequenceContent
        mod <- "Per_base_sequence_content"
        df <- getModule(x, mod)

        if (!length(df)) {
            p <- .emptyPlot("No Sequence Content Module Detected")
            if (usePlotly) p <- ggplotly(p, tooltip = "")
            return(p)
        }

        # Sort out any binned positions
        df$Base <- lapply(
            df$Base,
            function(x){
                rng <- as.integer(str_split(x, pattern = "-")[[1]])
                seq(min(rng), max(rng), by = 1L)
            }
        )
        df <- unnest(df, Base)

        plotType <- match.arg(plotType)
        if (missing(pwfCols)) pwfCols <- pwf
        stopifnot(is(plotTheme, "theme"))

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(x, labels, pattern = pattern, ...)
        labels <- labels[names(labels) %in% df$Filename]

        ## Define the bases as a vector for ease later in the function
        acgt <- c("T", "C", "A", "G")
        ## Axis labels
        xLab <- "Position in read (bp)"
        yLab <- "Percent (%)"

        ## Get the PASS/WARN/FAIL status
        status <- getSummary(x)
        status <- subset(status, Category == "Per base sequence content")
        status$Status <- factor(
            status$Status, levels = c("PASS", "WARN", "FAIL")
        )
        status <- droplevels(status)

        if (plotType == "heatmap") {

            ## Round to 2 digits to reduce the complexity of the colour
            ## palette
            df[acgt] <- lapply(df[acgt], round, digits = 2)
            maxBase <- max(vapply(acgt, function(x){max(df[[x]])}, numeric(1)))
            ## Set the colours, using opacity for G
            df$opacity <- 1 - df$G / maxBase
            df$RGB <- with(df, rgb(
                red = `T` * opacity / maxBase,
                green = A * opacity / maxBase,
                blue = C * opacity / maxBase)
            )

            basicStat <- getModule(x, "Basic_Statistics")
            basicStat <- basicStat[c("Filename", "Longest_sequence")]
            df <- dplyr::right_join(df, basicStat, by = "Filename")
            cols <- c("Filename", "Base", "RGB", "Longest_sequence")
            df <- df[c(cols, acgt)]

            ## Now define the order for a dendrogram if required
            key <- names(labels)
            df_long <- tidyr::pivot_longer(
                df, cols = all_of(acgt), names_to = "Nt", values_to = "Percent"
            )
            df_long$Base <- paste(df_long$Base, df_long$Nt)
            df_long <- df_long[c("Filename", "Base", "Percent")]
            clusterDend <- .makeDendro(df_long, "Filename", "Base", "Percent")
            dx <- ggdendro::dendro_data(clusterDend)
            if (dendrogram | cluster) key <- labels(clusterDend)
            if (!dendrogram) dx$segments <- dx$segments[0,]
            ## Now set everything as factors
            df$Filename <- factor(labels[df$Filename], levels = labels[key])
            status$Filename <- factor(
                labels[status$Filename], levels = labels[key]
            )

            ## Define the colours as named colours (name = RGB)
            tileCols <- unique(df$RGB)
            names(tileCols) <- unique(df$RGB)

            ## Define the tile locations
            df$y <- as.integer(df$Filename)
            df$ymax <- as.integer(df$Filename) + 0.5
            df$ymin <- df$ymax - 1
            df$xmax <- df$Base + 0.5
            df$xmin <- df$Base - 1
            ## Add percentage signs to ACGT for prettier labels
            df[acgt] <- lapply(
                df[acgt], scales::percent, accuracy = 0.1, scale = 1
            )
            yBreaks <- seq_along(levels(df$Filename))

            p <- ggplot(
                df,
                aes(
                    A = !!sym("A"), C = !!sym("C"),
                    G = !!sym("G"), `T` = !!sym("T"),
                    Filename = Filename, Position = Base, fill = !!sym("RGB")
                )
            ) +
                geom_rect(
                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                    linetype = 0
                ) +
                ggtitle("Per Base Sequence Content") +
                scale_fill_manual(values = tileCols) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(
                    expand = c(0, 0), breaks = yBreaks,
                    labels = levels(df$Filename),
                    position = "right"
                ) +
                labs(x = xLab, y = c()) +
                theme_bw() +
                theme(
                    legend.position = "none", panel.grid = element_blank(),
                    plot.title = element_text(
                        hjust = 0.5 * heat_w / (heat_w + 1)
                    )
                ) +
                plotTheme
            if (showPwf | dendrogram)
                p <- p +
                theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))

            if (!showPwf) status <- status[0,]

            p <- .prepHeatmap(
                p, status, dx$segments, usePlotly, heat_w, pwfCols
            )

        }

        if (plotType == "line") {

            df$Filename <- labels[df$Filename]
            df <- tidyr::gather(df, "Nt", "Percent", one_of(acgt))
            df$Nt <- factor(df$Nt, levels = acgt)
            df$Percent <- round(df$Percent, 2) / 100

            ## Add the pwf status for plotly plots
            status[["Filename"]] <- labels[status[["Filename"]]]
            df <- left_join(df, status, by = "Filename")
            rect_df <- group_by(df, Filename, Status)
            rect_df <- dplyr::summarise(
                rect_df, Start = 0, End = max(df$Base), .groups = "drop"
            )

            df$diff <- c(Inf, diff(df$Percent))
            df <- subset(df, diff != 0)

            ## Set colours for line plots & theme
            if (is.null(scaleColour)) {
                baseCols <- c(`T` = "red", G = "black", A = "green", C = "blue")
                scaleColour <- scale_colour_manual(values = baseCols)
            }
            stopifnot(is(scaleColour, "ScaleDiscrete"))
            stopifnot(scaleColour$aesthetics == "colour")

            p <- ggplot(df, aes(Base, Percent, colour = Nt))
            if (showPwf) p <- p + geom_rect(
                aes(xmin = 0, xmax = End, ymin = 0, ymax = 1, fill = Status),
                data = rect_df,
                alpha = 0.1, inherit.aes = FALSE
            )

            p <- p + geom_line(...) +
                facet_wrap(~Filename, ncol = nc) +
                scale_y_continuous(
                    limits = c(0, 1), expand = c(0, 0), labels = scales::percent
                ) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_fill_manual(
                    values = getColours(pwfCols)[levels(status$Status)]
                ) +
                scaleColour +
                labs(x = xLab, y = yLab, colour = "Base") +
                theme_bw() +
                plotTheme

            if (usePlotly) {
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                ttip <- c("y", "colour", "label", "fill")
                p <- suppressMessages(
                    suppressWarnings(plotly::ggplotly(p, tooltip = ttip))
                )
            }
        }
        if (plotType == "residuals") {

            df$Filename <- labels[df$Filename]
            ## Convert to long form
            df <- pivot_longer(
                data = df, cols = all_of(acgt), names_to = "Nt",
                values_to = "Percent"
            )
            ## Calculate the Residuals for each base/position
            df <- group_by(df, Base, Nt)
            df <- dplyr::mutate(df, Residuals = Percent - mean(Percent))
            df <- ungroup(df)
            df[["Residuals"]] <- round(df[["Residuals"]], 2)
            ## Find the duplicated positions as a result of binning & remove
            df <- dplyr::arrange(df, Filename, Nt, Base)
            df <- group_by(df, Filename, Nt)
            df <- dplyr::mutate(df, diff = c(0, diff(Percent)))
            df <- ungroup(df)
            df <- dplyr::filter(df, diff != 0 | Base == 1)
            df[["Deviation"]] <- percent(df[["Residuals"]]/100, accuracy = 0.1)
            ## Add the pwf status for plotly plots
            status[["Filename"]] <- labels[status[["Filename"]]]
            df <- left_join(df, status, by = "Filename")

            ## Set colours for line plots & theme
            if (is.null(scaleColour)) {
                scaleColour <- scale_colour_brewer(palette = "Paired")
            }
            stopifnot(is(scaleColour, "ScaleDiscrete"))
            stopifnot(scaleColour$aesthetics == "colour")

            Deviation <- c()
            p <- ggplot(
                df,
                aes(
                    Base, Residuals, colour = Filename, label = Deviation,
                    status = Status
                )
            ) +
                geom_line(aes(group = Filename), ...) +
                facet_wrap(~Nt) +
                scale_y_continuous(labels = .addPercent) +
                scale_x_continuous(expand = c(0, 0)) +
                scaleColour +
                labs(x = xLab) +
                theme_bw() +
                plotTheme

            if (usePlotly) {
                ttip <- c("x", "colour", "label", "status")
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                p <- suppressMessages(
                    suppressWarnings(plotly::ggplotly(p, tooltip = ttip))
                )
            }

        }

        p
    }
)
#' @importFrom tidyr unnest pivot_longer
#' @importFrom rlang sym !!
#' @importFrom scales label_percent
#' @rdname plotSeqContent-methods
#' @export
setMethod(
    "plotSeqContent", signature = "FastpData",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        module = c("Before_filtering", "After_filtering"),
        reads = c("read1", "read2"), readsBy = c("facet", "linetype"),
        moduleBy = c("facet", "linetype"),
        bases = c("A", "T", "C", "G", "N", "GC"), scaleColour = NULL,
        scaleLine = NULL, plotlyLegend = FALSE, plotTheme = theme_get(),
        expand.x = 0.02, expand.y = c(0, 0.05), ...
    ) {

        module <- match.arg(module, several.ok = TRUE)
        nMod <- length(module)
        reads <- match.arg(reads, several.ok = TRUE)
        nReads <- length(reads)
        data <- lapply(module, function(i) getModule(x, i)[reads])
        names(data) <- module
        data <- lapply(data, dplyr::bind_rows)
        df <- bind_rows(data, .id = "Module")
        df <- df[c("Filename", "fqName", "Module", "content_curves")]
        df <- unnest(df, !!sym("content_curves"))
        bases <- match.arg(bases, colnames(df), several.ok = TRUE)
        df <- pivot_longer(
            df, all_of(bases), names_to = "Base", values_to = "Frequency"
        )
        df$Base <- factor(df$Base, levels = bases)
        df$Module <- factor(df$Module, levels = module)

        ## Sort out plotting args
        readsBy <- match.arg(readsBy)
        moduleBy <- match.arg(moduleBy)
        if (readsBy == moduleBy & readsBy != "facet") stop(
            "Cannot set the same plotting parameter to both reads and module"
        )
        lt <- NULL
        if (readsBy == "linetype" & nReads > 1) lt <- sym("fqName")
        if (moduleBy == "linetype" & nMod > 1) lt <- sym("Module")
        fm <- dplyr::case_when(
            readsBy == "facet" & moduleBy == "facet" ~ "Module ~ fqName",
            readsBy != "facet" & moduleBy == "facet" ~ ". ~ Module",
            readsBy == "facet" & moduleBy != "facet" ~ ". ~ fqName",
            TRUE ~ ". ~ ."
        )
        fm <- as.formula(fm)

        if (is.null(scaleColour)) {
            ## Best guess based on those in a fastp report
            basecols <- c("#807C58", "#601490", "green", "blue", "red","grey20")
            names(basecols) <- c("A", "T", "C", "G", "N", "GC")
            scaleColour <- scale_colour_manual(values = basecols[bases])
        }
        stopifnot(is(scaleColour, "ScaleDiscrete"))
        stopifnot(scaleColour$aesthetics == "colour")
        if (is.null(scaleLine)) scaleLine <- scale_linetype_discrete()
        stopifnot(is(scaleLine, "ScaleDiscrete"))
        stopifnot(scaleLine$aesthetics == "linetype")
        stopifnot(is(plotTheme, "theme"))

        ## Sort out labels for nicer plotting
        fqLabels <- .makeLabels(df, pattern = pattern, col = "fqName")
        df$fqName <- factor(fqLabels[df$fqName], levels = fqLabels)
        labels <- .makeLabels(df, labels, pattern)
        df$Filename <- factor(labels[df$Filename], levels = labels)
        main <- unique(labels)

        ## Final tweaks for better plotly category labels
        names(df) <- gsub("position", "Position", names(df))
        df$Frequency <- round(100 * df$Frequency, 2)
        df$Filename <- df$fqName
        p <- ggplot(
            df,
            aes(
                Position, Frequency, colour = Base, linetype = {{ lt }},
                name = Filename
            )
        ) +
            geom_line(...) +
            facet_grid(fm) +
            ggtitle(main) +
            labs(x = "Position in Read") +
            scale_x_continuous(expand = expansion(rep_len(expand.x, 2))) +
            scale_y_continuous(
                labels = label_percent(scale = 1),
                limits = c(0, max(df$Frequency)),
                expand = expansion(rep_len(expand.y, 2))
            ) +
            scaleColour + scaleLine +
            theme_bw() + plotTheme

        if (usePlotly) {
            if (!plotlyLegend) p <- p + theme(legend.position = "none")
            tt <- c("Filename", "Position", "Frequency", "Base")
            p <- plotly::ggplotly(p, tooltip = tt)
        }
        p
    }
)
#' @importFrom dplyr bind_rows distinct summarise case_when left_join group_by
#' @importFrom tidyr unnest pivot_longer
#' @importFrom tidyselect all_of
#' @rdname plotSeqContent-methods
#' @export
setMethod(
    "plotSeqContent", "FastpDataList",
    function(
        x, usePlotly = FALSE, labels, pattern = ".(fast|fq|bam).*",
        module = c("Before_filtering", "After_filtering"),
        moduleBy = c("facet", "linetype"),
        reads = c("read1", "read2"), readsBy = c("facet", "linetype"),
        bases = c("A", "T", "C", "G", "N", "GC"), showPwf = FALSE, pwfCols,
        warn = 10, fail = 20, plotType = c("heatmap", "line", "residuals"),
        plotlyLegend = FALSE, scaleColour = NULL, scaleLine = NULL,
        plotTheme = theme_get(),
        cluster = FALSE, dendrogram = FALSE, heat_w = 8,
        expand.x = c(0.01), expand.y = c(0, 0.05), nc = 2, ...
    ){

        ## Check args
        # We can't plot B4/After if we cluster
        mod <- match.arg(module, several.ok = TRUE)
        nMod <- length(mod)
        reads <- match.arg(reads, several.ok = TRUE)
        nReads <- length(reads)

        ## Setup the data
        data <- lapply(mod, function(m) getModule(x, m)[reads])
        names(data) <- mod
        data <- lapply(data, bind_rows, .id = "reads")
        df <- bind_rows(data, .id = "Module")
        df[["Module"]] <- factor(df[["Module"]], levels = mod)
        df <- df[c("Filename", "fqName", "Module", "reads", "content_curves")]
        df <- unnest(df, !!sym("content_curves"))
        bases <- match.arg(bases, colnames(df), several.ok = TRUE)

        if (!length(df)) {
            p <- .emptyPlot("No Sequence Content Module Detected")
            if (usePlotly) p <- plotly::ggplotly(p, tooltip = "")
            return(p)
        }

        ## Set any Pwf Status
        status_df <- summarise(
            df,
            AT = max(abs(!!sym("A") - !!sym("T"))),
            GC = max(abs(!!sym("G") - !!sym("C"))), .by = Filename
        )
        status_df$Status <- case_when(
            status_df$AT > fail / 100 | status_df$GC > fail / 100 ~ "FAIL",
            status_df$AT > warn / 100 | status_df$GC > warn / 100 ~ "WARN",
            TRUE ~ "PASS"
        )
        status_df$Status <- factor(
            status_df$Status, levels = slotNames("PwfCols")[seq_len(3)]
        )
        status_df <- status_df[c("Filename", "Status")]
        plotType <- match.arg(plotType)
        if (missing(pwfCols)) pwfCols <- pwf
        stopifnot(is(plotTheme, "theme"))

        ## Drop the suffix, or check the alternate labels
        fqLabels <- .makeLabels(df, labels, pattern = pattern, col = "fqName")
        df$fqName <- factor(fqLabels[df$fqName], levels = fqLabels)
        labels <- .makeLabels(df, labels, pattern = pattern)
        key <- names(labels)

        ## Layout for line plots
        readsBy <- match.arg(readsBy)
        moduleBy <- match.arg(moduleBy)
        if (
            all(
                readsBy == moduleBy, plotType != "heatmap", nMod > 1,
                nReads > 1
            )
        ) stop(
            "Cannot set reads and module to the same parameter for line plots"
        )
        linetype <- NULL
        if (readsBy == "linetype" & nReads > 1) linetype <- sym("reads")
        if (moduleBy == "linetype" & nMod > 1) linetype <- sym("Module")
        if (is.null(scaleLine)) scaleLine <- scale_linetype_discrete()
        stopifnot(is(scaleLine, "ScaleDiscrete"))
        stopifnot(scaleLine$aesthetics == "linetype")

        ## Axis labels
        xLab <- "Position in read (bp)"
        yLab <- "Percent (%)"
        main <- paste("Sequence Content:", paste(reads, collapse = " & "))
        main <- stringr::str_to_title(main)

        if (plotType == "heatmap") {

            n_facets <- length(mod) * length(reads)
            if ((dendrogram | cluster) & n_facets > 2) {
                message(
                    "Cannot cluster when plotting both modules and both sets ",
                    "of reads. Setting cluster and dendrogram to FALSE"
                )
                cluster <- dendrogram <- FALSE
            }

            if (any(c("N", "GC") %in% bases))
                message("N/GC bases will be ignored when preparing a heatmap")
            bases <- c("A", "C", "G", "T")
            df <- dplyr::select(
                df, Filename, !!sym("fqName"), !!sym("Module"),
                !!sym("reads"), all_of(bases), !!sym("position")
            )
            df[bases] <- lapply(df[bases], function(x) round(100 * x, 2))

            ## Round to 2 digits to reduce the complexity of the colour palette
            maxFreq <- max(unlist(df[bases]))
            ## Set the colours, using opacity for G
            df$opacity <- 1 - df$G / maxFreq
            df$RGB <- rgb(
                red = df[["T"]] * df$opacity / maxFreq,
                green = df[["A"]] * df$opacity / maxFreq,
                blue = df[["C"]] * df$opacity / maxFreq
            )
            names(df) <- gsub("position", "Position", names(df))

            ## Set up the dendrogram & labels
            n <- length(labels)
            if (n == 1 & (cluster | dendrogram))
                message(
                    "Cannot cluster one file. Ignoring cluster and dendgrogram"
                )
            ## Now define the order for a dendrogram if required
            if (n > 1 & (cluster | dendrogram)) {
                df_long <- tidyr::pivot_longer(
                    df, cols = all_of(bases), names_to = "Nt",
                    values_to = "Percent"
                )
                df_long$Position <- paste(
                    df_long$Position, df_long$Nt, df_long$reads, df_long$Module
                )
                df_long <- df_long[c("Filename", "Position", "Percent")]
                clusterDend <- .makeDendro(
                    df_long, "Filename", "Position", "Percent"
                )
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
            df$Filename <- factor(labels[df$Filename], levels = labels[key])
            status_df$Filename <- factor(
                labels[status_df$Filename], levels = labels[key]
            )

            ## Define the colours as named colours (name = RGB)
            tileCols <- unique(df$RGB)
            names(tileCols) <- unique(df$RGB)

            ## Define the tile locations
            df$y <- as.integer(df$Filename)
            df$ymax <- as.integer(df$Filename) + 0.5
            df$ymin <- df$ymax - 1
            df$xmax <- df$Position + 0.5
            df$xmin <- df$Position - 1
            ## Add percentage signs to ACGT for prettier labels
            df[bases] <- lapply(
                df[bases], scales::percent, accuracy = 0.1, scale = 1
            )
            yBreaks <- seq_along(levels(df$Filename))

            fm <- case_when(
                n_facets == 1 ~ ". ~ .",
                n_facets == 2 & length(mod) == 2 ~ ". ~ Module",
                n_facets == 2 & length(reads) == 2 ~ ". ~ reads",
                TRUE ~ "Module ~ reads"
            )
            fm <- as.formula(fm)

            p <- ggplot(
                df,
                aes(
                    A = !!sym("A"), C = !!sym("C"),
                    G = !!sym("G"), `T` = !!sym("T"),
                    fqName = !!sym("fqName"), Position = Position,
                    fill = !!sym("RGB")
                )
            ) +
                geom_rect(
                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                    linetype = 0
                ) +
                facet_grid(fm, switch = "y") +
                ggtitle(main) +
                scale_fill_manual(values = tileCols) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(
                    expand = c(0, 0), breaks = yBreaks,
                    labels = levels(df$Filename),
                    position = "right"
                ) +
                labs(x = xLab, y = c()) +
                theme_bw() +
                theme(
                    legend.position = "none", panel.grid = element_blank(),
                    plot.title = element_text(
                        hjust = 0.5 * heat_w / (heat_w + 1)
                    )
                ) +
                plotTheme
            if (showPwf | dendrogram)
                p <- p + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"))
            tt <- c(bases, "fqName", "Position")
            if (!showPwf) status_df <- status_df[0, ]
            p <- .prepHeatmap(
                p, status_df, dx$segments, usePlotly, heat_w, pwfCols, tt
            )

        }

        if (plotType == "line") {

            df$Filename <- factor(labels[df$Filename], labels)
            df <- pivot_longer(
                df, all_of(bases), names_to = "Base", values_to = "Frequency"
            )
            df$Base <- factor(df$Base, levels = bases)
            df$Frequency <- round(100 * df$Frequency, 2)

            ## Add the pwf status for plotly plots
            status_df$Filename <- factor(
                labels[status_df$Filename], levels = labels[key]
            )
            df <- left_join(df, status_df, by = "Filename")
            rect_df <- distinct(df, Filename, Status)
            expand.x <- rep_len(expand.x, 2)
            x_lim <- c(0, max(df$position))
            rng <- diff(range(x_lim))
            x_lim <- x_lim + c(-1, 1) * expand.x * rng
            rect_df$Start <- x_lim[[1]]
            rect_df$End <- x_lim[[2]]

            ## Set colours for line plots & theme
            if (is.null(scaleColour)) {
                baseCols <- c(
                    "#807C58", "#601490", "green", "blue", "red", "grey20"
                )
                names(baseCols) <- c("A", "T", "C", "G", "N", "GC")
                scaleColour <- scale_colour_manual(values = baseCols[bases])
            }
            stopifnot(is(scaleColour, "ScaleDiscrete"))
            stopifnot(scaleColour$aesthetics == "colour")

            names(df) <- gsub("position", "Position", names(df))
            expand.y <- rep_len(expand.y, 2)
            y_lim <- c(0, max(df$Frequency))
            rng <- diff(range(y_lim))
            y_lim <- y_lim + c(-1, 1) * expand.y * rng

            ## These settings are specific to plotType == "line" so should stay here
            fm <- dplyr::case_when(
                readsBy == "facet" & nReads > 1 ~ "Filename ~ reads",
                moduleBy == "facet" & nMod > 1 ~ "Filename ~ Module",
                TRUE ~ "Filename ~ ."
            )
            facets <- facet_grid(as.formula(fm))

            p <- ggplot(
                df,
                aes(
                    Position, Frequency, colour = Base,
                    linetype = {{ linetype }}, reads = !!sym("fqName")
                )
            )
            ## Add rectangles if requested
            if (showPwf) p <- p + geom_rect(
                aes(
                    xmin = Start, xmax = End, ymin = 0, ymax = max(y_lim),
                    fill = Status
                ),
                data = rect_df, alpha = 0.1, inherit.aes = FALSE
            )
            ## The rest of the plot
            p <- p + geom_line(...) + facets +
                scale_y_continuous(
                    limits = y_lim, expand = rep_len(0, 4),
                    labels = scales::label_percent(scale = 1)
                ) +
                scale_x_continuous(limits = x_lim, expand = rep_len(0, 4)) +
                scale_fill_manual(
                    values = getColours(pwfCols)[levels(status_df$Status)]
                ) +
                scaleColour + scaleLine +
                ggtitle(main) +
                labs(x = xLab, y = yLab, colour = "Base") +
                theme_bw() + plotTheme

            if (usePlotly) {
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                ttip <- c("y", "colour", "label", "fill", "fqName")
                p <- suppressMessages(
                    suppressWarnings(plotly::ggplotly(p, tooltip = ttip))
                )
            }
        }

        if (plotType == "residuals") {

            df$Filename <- labels[df$Filename]
            ## Convert to long form
            df <- pivot_longer(
                df, cols = all_of(bases), names_to = "Base",
                values_to = "Percent"
            )
            names(df) <- gsub("position", "Position", names(df))
            ## Calculate the Residuals for each base/position
            df <- group_by(df, Position, Base, reads)
            df <- dplyr::mutate(df, Residuals = Percent - mean(Percent))
            df <- ungroup(df)
            df[["Residuals"]] <- round(100 * df[["Residuals"]], 2)
            df[["Deviation"]] <- percent(df[["Residuals"]]/100, accuracy = 0.1)

            ## Add the pwf status for plotly plots
            status_df[["Filename"]] <- labels[status_df[["Filename"]]]
            df <- left_join(df, status_df, by = "Filename")

            ## Set colours for line plots & theme
            if (is.null(scaleColour)) {
                n <- length(labels)
                if (n <= 9) {
                    scaleColour <- scale_colour_brewer(palette = "Set1")
                } else {
                    cols <- hcl.colors(n, "Dark 2")
                    names(cols) <- labels
                    scaleColour <- scale_colour_manual(values = cols)
                }
            }
            stopifnot(is(scaleColour, "ScaleDiscrete"))
            stopifnot(scaleColour$aesthetics == "colour")

            ## These settings are specific to plotType == "line" so should stay here
            fm <- dplyr::case_when(
                readsBy == "facet" & nReads > 1~ "Base ~ reads",
                moduleBy == "facet" & nMod > 1 ~ "Base ~ Module",
                TRUE ~ "Base ~ ."
            )
            facets <- facet_grid(as.formula(fm))

            p <- ggplot(
                df,
                aes(
                    Position, Residuals, colour = Filename,
                    label = !!sym("Deviation"), status = Status,
                    fqName = !!sym("fqName"), linetype = {{ linetype }}
                )
            ) +
                geom_line(...) + facets +
                scale_y_continuous(labels = label_percent(scale = 1)) +
                scale_x_continuous(expand = expansion(rep_len(expand.x, 2))) +
                scaleColour + scaleLine +
                ggtitle(main) +
                labs(x = xLab) +
                theme_bw() + plotTheme

            if (usePlotly) {
                ttip <- c("x", "colour", "label", "linetype", "fqName")
                if (showPwf) ttip <- c(ttip, "status")
                if (!plotlyLegend) p <- p + theme(legend.position = "none")
                p <- suppressMessages(
                    suppressWarnings(plotly::ggplotly(p, tooltip = ttip))
                )
            }

        }

        p
    }
)
