#' @title Estimate a GC Content Distribution From Sequences
#'
#' @description Generate a GC content distribution from sequences for a given
#' read length and fragment length
#'
#' @param x `DNAStringSet` or path to a fasta file
#' @param n The number of reads to sample
#' @param rl Read Lengths to sample
#' @param fl The mean of the fragment lengths sequenced
#' @param fragSd The standard deviation of the fragment lengths being sequenced
#' @param bins The number of bins to estimate
#' @param ... Not used
#'
#' @details
#' The function takes the supplied object and returns the theoretical GC
#' content distribution. Using a fixed read length essentially leads to a
#' discrete distribution so the bins argument is used to define the number of
#' bins returned. This defaults to 101 for 0 to 100% inclusive.
#'
#' The returned values are obtained by interpolating the values obtained during
#' sampling. This avoids returned distributions with gaps and jumps as would be
#' obtained setting readLengths at values not in multiples of 100.
#'
#' Based heavily on https://github.com/mikelove/fastqcTheoreticalGC
#'
#' @return A `tibble` with two columns: `GC_Content` and `Freq`
#' denoting the proportion of GC and frequency of occurence reqpectively
#'
#' @examples
#' faDir <- system.file("extdata", package = "ngsReports")
#' faFile <- list.files(faDir, pattern = "fasta", full.names = TRUE)
#' df <- estGcDistn(faFile, n = 200)
#'
#' @docType methods
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet subseq letterFrequency
#' @importFrom BiocGenerics width
#' @importFrom stats rnorm runif lm.fit
#'
#' @export
#' @rdname estGcDistn
setGeneric(
    "estGcDistn",
    function(
        x, n = 1e6, rl = 100, fl = 200, fragSd = 30, bins = 101, ...
    ) standardGeneric("estGcDistn")
)
#' @export
#' @rdname estGcDistn
#' @aliases estGcDistn
setMethod("estGcDistn", "ANY", function(x, ...){
    cl <- class(x)
    message("No method defined for objects of class ", cl)
})
#' @export
#' @rdname estGcDistn
#' @aliases estGcDistn
setMethod(
    "estGcDistn", "character",
    function(x, n = 1e6, rl = 100, fl = 200, fragSd = 30, bins = 101, ...) {

        stopifnot(file.exists(x))
        errMsg <- "Supplied File not in FASTA format"
        x <- tryCatch(
            readDNAStringSet(x, format = "fasta"),
            error = function(e) {stop(errMsg)}
        )
        estGcDistn(x, n, rl, fl, fragSd, bins)
    })
#' @export
#' @rdname estGcDistn
#' @aliases estGcDistn
setMethod("estGcDistn", "DNAStringSet", function(
        x, n = 1e6, rl = 100, fl = 200, fragSd = 30, bins = 101, ...) {

    if (!requireNamespace("truncnorm", quietly = TRUE)) {
        stop(
            "The package truncnorm is required for this function.\n",
            "Please install it using BiocManager::install('truncnorm')"
        )
    }

    ## Check the arguments & convert to integers
    n <- tryCatch(as.integer(n))
    rl <- tryCatch(as.integer(rl))
    fl <- tryCatch(as.integer(fl))
    bins <- tryCatch(as.integer(bins))
    stopifnot(rl <= fl)
    stopifnot(is.numeric(fragSd))
    stopifnot(fragSd >= 0)

    ## Randomly select sequences for the fragments
    idx <- sample(seq_along(x), n, replace = TRUE)
    molecules <- x[idx]

    ## Sample a set of fragments with variable length using the truncated norm
    ## Upper & lower limits are set as 1 & sequence length respectively
    fragSizes <- truncnorm::rtruncnorm(n, 1, width(molecules), fl, fragSd)
    fragSizes <- floor(fragSizes)

    ## sample start positions uniformly
    starts <- runif(n, 1, width(molecules) - fragSizes)
    starts <- floor(starts)

    ## Set the end points
    ends <- starts + fragSizes - 1
    ends[ends > width(molecules)] <- width(molecules)[ends > width(molecules)]

    ## Generate the fragments
    frags <- subseq(molecules, start = starts, end = ends)

    ## Now sample down to reads.
    ## Where frags are < readLength, sample the shorter value
    rdEnds <- width(frags)
    rdEnds[rdEnds > rl] <- rl
    reads <- subseq(frags, start = 1, end = rdEnds)

    ## Find the GC content for each read
    gc <- as.vector(letterFrequency(reads, "GC"))
    total <- as.vector(letterFrequency(reads, "ACGT"))
    gc.content <- (gc / total)[total > 0]

    ## Form into bins based on RL
    breaks <- seq(0, rl) / rl
    gcCut <- table(cut(gc.content, breaks = breaks, include.lowest = FALSE))
    freq <- c(sum(gc.content == 0), gcCut) / length(frags)

    ## Now we'll interpolate using sets of three points to fit regression lines
    ## This deals with the issues trying to obtain a continuous distribution
    ## when dividing by a discrete denominator
    ## Unfortunately this is necessary as FastQC only provides values for
    ## percentages from 0 to 100 in integer steps
    splits <- seq(1, length(breaks) - 2)
    lmFits <- lapply(splits, function(x){
        vals <- seq(x, x + 2)
        fit <- lm.fit(x = cbind(1, breaks[vals]), y = freq[vals])
        fit$coefficients
    })

    ## Setup the values to return
    props <- seq(0, by = 1, length.out = bins) / (bins - 1)
    df <- tibble(
        GC_Content = props * 100,
        Freq = vapply(props, function(x){
            interval <- findInterval(x, breaks[splits])
            max(lmFits[[interval]]["x2"]*x + lmFits[[interval]]["x1"], 0)
        }, numeric(1)))

    df
}
)
