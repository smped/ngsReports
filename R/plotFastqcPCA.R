#' @title Draw a PCA plot for Fast QC modules
#'
#' @description Draw a PCA plot for Fast QC modules across multiple samples
#' \lifecycle{experimental}
#'
#' @details
#' This carries out PCA on a single FastQC module and plots the
#' output using either ggplot or plotly. Current modules for PCA are
#' Per_base_sequence_quality, Per_sequence_quality_scores,
#' Per_sequence_GC_content, Per_base_sequence_content, and
#' Sequence_Length_Distribution.
#'
#' If a factor is provided in the groups argument, this will be applied to the
#' plotting colours and ellipses will be drawn using these groups.
#' Only the labels will be plotted using `geom_text()`
#'
#' @param x Can be a `FastqcDataList` or `character` vector of file paths
#' @param module `character` vector containing
#'  the desired FastQC module (eg. c("Per_base_sequence_quality",
#'  "Per_base_sequence_content"))
#' @param usePlotly `logical`. Output as ggplot2 (default) or plotly
#' object.
#' @param labels An optional named vector of labels for the file names.
#' All file names must be present in the names of the vector.
#' @param pattern Regex to remove from the end of any filenames
#' @param sz The size of the text labels
#' @param pc The two components to be plotted
#' @param groups Optional factor of the same length as x. If provided,
#' the plot will be coloured using this factor as the defined groups.
#' Ellipses will also be added to the final plot.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @docType methods
#'
#' @importFrom stats prcomp
#' @importFrom scales percent
#' @importFrom tidyr pivot_wider
#' @importFrom zoo na.locf
#' @importFrom rlang !! sym
#' @import ggplot2
#' @import tibble
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#' grp <- as.factor(gsub(".+(R[12]).*", "\\1", fqName(fdl)))
#' plotFastqcPCA(fdl, module = "Per_sequence_GC_content", groups = grp)
#'
#'
#' @name plotFastqcPCA
#' @rdname plotFastqcPCA-methods
#' @export
setGeneric(
    "plotFastqcPCA",
    function(
        x, module = "Per_sequence_GC_content", usePlotly = FALSE, labels,
        pattern = ".(fast|fq|bam).*", sz = 4, groups, ...
    ) standardGeneric("plotFastqcPCA")
)
#' @rdname plotFastqcPCA-methods
#' @export
setMethod(
    "plotFastqcPCA", signature = "ANY",
    function(
        x, module = "Per_sequence_GC_content", usePlotly = FALSE, labels,
        pattern = ".(fast|fq|bam).*", sz = 4, groups, ...
    ){.errNotImp(x)}
)
#' @rdname plotFastqcPCA-methods
#' @export
setMethod(
    "plotFastqcPCA", signature = "character",
    function(
        x, module = "Per_sequence_GC_content", usePlotly = FALSE, labels,
        pattern = ".(fast|fq|bam).*", sz = 4, groups, ...
    ){
        x <- FastqcDataList(x)
        if (length(x) == 1) x <- FastqcData(x[[1]]) # This will raise an error
        plotFastqcPCA(x, module, usePlotly, labels, pattern, sz, groups, ... )
    }
)
#' @rdname plotFastqcPCA-methods
#' @export
setMethod(
    "plotFastqcPCA", signature = "FastqcDataList",
    function(
        x, module = "Per_sequence_GC_content", usePlotly = FALSE, labels,
        pattern = ".(fast|fq|bam).*", sz = 4, groups, pc = c(1, 2), ...
    ){

        availMods <- c(
            "Per_base_sequence_quality", "Per_sequence_quality_scores",
            "Per_sequence_GC_content", "Per_base_sequence_content",
            "Sequence_Length_Distribution"
        )
        module <- match.arg(module, availMods)
        ## Ensure two components are selected as integers
        pc <- as.integer(pc)
        stopifnot(length(unique(pc)) > 1)
        stopifnot(!any(is.na(pc)))
        pc <- paste0("PC", unique(pc)[seq_len(2)])


        ## Define the data preparation function
        pFun <- paste0(".generate", module, "PCA")
        args <- list(x = x)
        mat <- do.call(pFun, args)

        ## The PCA
        pca <- prcomp(mat, center = TRUE, scale. = TRUE)
        vars <- summary(pca)$importance["Proportion of Variance", pc]
        df <- as.data.frame(pca$x[, pc])
        df <- rownames_to_column(df, "Filename")

        ## Check any groups
        plot_aes <- aes(!!sym(pc[[1]]), !!sym(pc[[2]]))
        showGroups <- !missing(groups)
        if (showGroups) {
            stopifnot(length(groups) == length(x))
            stopifnot(is.factor(groups))
            df$group <- groups
            plot_aes <- aes(!!sym(pc[[1]]), !!sym(pc[[2]]), colour = group)
        }
        df$Label <- .makeLabels(df, labels, pattern)

        p <- ggplot(data = df, mapping = plot_aes) +
            geom_text(aes(label = Label), show.legend = FALSE, size = sz) +
            theme_bw() +
            labs(
                x = paste0(pc[[1]], " (", percent(vars[[1]]), ")"),
                y = paste0(pc[[2]], " (", percent(vars[[2]]), ")"),
                colour = "Group"
            )
        p <- .updateThemeFromDots(p, ...length())
        if (showGroups) p <- p + stat_ellipse()
        if (usePlotly) p <- plotly::ggplotly(p)
        p
    }

)


.generatePer_base_sequence_qualityPCA <- function(x){

    df <- getModule(x, "Per_base_sequence_quality")
    df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))

    ## Adjust the data for files with varying read lengths
    ## This will fill NA values with the previous values
    df <- lapply(
        split(df, f = df$Filename),
        function(y){
            Longest_sequence <-
                gsub(".*-([0-9]*)", "\\1", as.character(y$Base))
            Longest_sequence <- max(as.integer(Longest_sequence))
            dfFill <- data.frame(Start = seq_len(Longest_sequence))
            y <- dplyr::right_join(y, dfFill, by = "Start")
            na.locf(y)
        }
    )

    df <- dplyr::bind_rows(df)[c("Filename", "Start", "Mean")]
    df <- pivot_wider(
        data = df, names_from = "Start", values_from = "Mean", values_fill = 0
    )
    df <- as.data.frame(df)
    df <- column_to_rownames(df, "Filename")
    df
}

.generatePer_sequence_quality_scoresPCA <- function(x){

    df <- getModule(x, "Per_sequence_quality_scores")
    df <- pivot_wider(
        df, names_from = "Quality", values_from = "Count", values_fill = 0
    )
    df <- as.data.frame(df)
    df <- column_to_rownames(df, "Filename")
    df
}


.generatePer_sequence_GC_contentPCA <- function(x){

    df <- getModule(x, "Per_sequence_GC_content")
    df <- pivot_wider(
        df, names_from = "GC_Content", values_from = "Count", values_fill = 0
    )
    df <- as.data.frame(df)
    df <- column_to_rownames(df, "Filename")
    df
}


.generatePer_base_sequence_contentPCA <- function(x){

    df <- getModule(x, "Per_base_sequence_content")
    df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))

    ## Adjust the data for files with varying read lengths
    ## This will fill NA values with the previous values
    df <- lapply(split(df, f = df$Filename), function(y){
        Longest_sequence <-
            gsub(".*-([0-9]*)", "\\1", as.character(y$Base))
        Longest_sequence <- max(as.integer(Longest_sequence))
        dfFill <- data.frame(Start = seq_len(Longest_sequence))
        y <- dplyr::right_join(y, dfFill, by = "Start")
        na.locf(y)
    })
    df <- dplyr::bind_rows(df)[c("Filename", "Start", "G", "A", "T", "C")]
    df <- pivot_wider(
        df, values_from = c("A", "C", "G", "T"), names_from = Start,
        values_fill = 0
    )

    df <- column_to_rownames(df, "Filename")
    df
}

.generateSequence_Length_DistributionPCA <- function(x){

    df <- getModule(x, "Sequence_Length_Distribution")
    df <- lapply(split(df, f = df$Filename), function(y){
        Longest_sequence <-
            gsub(".*-([0-9]*)", "\\1", as.character(y$Length))
        Longest_sequence <- max(as.integer(Longest_sequence))
        dfFill <- data.frame(Lower = seq_len(Longest_sequence))
        y <- dplyr::right_join(y, dfFill, by = "Lower")
        na.locf(y)
    })
    df <- dplyr::bind_rows(df)[c("Filename", "Lower", "Count")]
    df <- pivot_wider(
        df, names_from = "Lower", values_from = "Count", values_fill =  0
    )
    df <- column_to_rownames(df, "Filename")
    df
}
