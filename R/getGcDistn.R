#' @title Generate a GC Content Distribution From Sequences
#'
#' @description Generate a GC Content Distribution From Sequences
#'
#' @param x \code{DNAStringSet} or path to a fasta file
#' @param n The number of reads to sample
#' @param readLength Read Lengths to sample
#' @param fragLength The mean of the fragment lengths sequenced
#' @param fragSd The standard deviation of the fragment lengths being sequenced
#' @param bins The number of bins to estimate
#'
#' @details
#' The function takes the supplied object and returns the theoretical GC
#' content distribution. Using a fixed read length essentially leads to a
#' discrete distribution so the bins argument is used to define the number of
#' bins returned. This defaults to 101 for 0 to 100% inclusive.
#'
#' The returned values are obtained by interpolating the values obtained during
#' sampling. This avoids returned distributions with gaps and jumps as would be
#' obtained setting readLengths at values < 100.
#'
#' Based heavily on https://github.com/mikelove/fastqcTheoreticalGC
#'
#' @return A \code{tibble} with two columns: \code{GC_Content} and \code{Freq}
#' denoting the proportion of GC and frequency of occurence reqpectively
#'
#' @examples
#' faDir <- system.file("extdata", package = "ngsReports")
#' faFile <- list.files(faDir, pattern = "fasta", full.names = TRUE)
#' gen_df <- getGcDistn(faFile, n = 200)
#'
#' @importFrom stats rnorm runif lm.fit
#' @importFrom XVector subseq
#' @importFrom truncnorm rtruncnorm
#' @importClassesFrom Biostrings DNAStringSet
#'
#' @export
#' @rdname getGcDistn-methods
setGeneric("getGcDistn", function(
    x, n = 1e6, readLength = 100, fragLength = 200, fragSd = 30, bins = 101
){
    standardGeneric("getGcDistn")
})
#' @rdname getGcDistn-methods
#' @aliases getGcDistn,character
#' @export
setMethod("getGcDistn", "character", function(
    x, n = 1e6, readLength = 100, fragLength = 200, fragSd = 30,
    bins = 101){

    stopifnot(file.exists(x))
    x <- tryCatch(readDNAStringSet(x))
    getGcDistn(x, n, readLength, fragLength, fragSd, bins)

})
#' @rdname getGcDistn-methods
#' @aliases getGcDistn,DNAStringSet
#' @export
setMethod("getGcDistn", "DNAStringSet", function(
    x, n = 1e6, readLength = 100, fragLength = 200, fragSd = 30,
    bins = 101){

    ## Check the arguments & convert to integers
    n <- tryCatch(as.integer(n))
    readLength <- tryCatch(as.integer(readLength))
    fragLength <- tryCatch(as.integer(fragLength))
    bins <- tryCatch(as.integer(bins))
    stopifnot(readLength <= fragLength)
    stopifnot(is.numeric(fragSd))
    stopifnot(fragSd > 0)

    ## Randomly select sequences for the fragments
    idx <- sample(seq_along(x), n, replace = TRUE)
    molecules <- x[idx]
    ## Sample a set of fragments with variable length uing the truncated norm
    fragSizes <- rtruncnorm(n, 1, width(molecules), fragLength, fragSd)
    fragSizes <- floor(fragSizes)
    ## sample start positions uniformly
    starts <- runif(n, 1, width(molecules) - fragSizes)
    starts <- floor(starts)
    ## Find any molecules shorter than the sample fragment length
    keep <- (width(molecules) - starts) > fragSizes
    if (any(!keep)) message(paste(sum(!keep), "fragments will be skipped"))
    molecules <- molecules[keep]
    fragSizes <- fragSizes[keep]
    starts <- starts[keep]
    ## Generate the reads
    frags <- subseq(molecules, start = starts, width = fragSizes)

    ## Now sample reads. Where frags are < readLength sample the shorter value
    ends <- rep(readLength, sum(keep))
    ## Where RL > the sampled FL, keep the bases up to the sampled FL
    shortFrags <- ends > width(frags)
    ends[shortFrags] <- width(frags)[shortFrags]
    reads <- subseq(frags, start = 1, end = ends)

    ## Find the GC content for each read
    gc <- as.vector(Biostrings::letterFrequency(reads, "GC"))
    total <- as.vector(Biostrings::letterFrequency(reads, "ACGT"))
    gc.content <- (gc / total)[total > 0]

    ## Form into bins based on RL
    breaks <- 0:readLength / readLength
    gcCut <- table(cut(gc.content, breaks = breaks, include.lowest = FALSE))
    freq <- c(sum(gc.content == 0), gcCut) / sum(keep)
    ## Now we'll interpolate using sets of three points to fit regression lines
    ## This deals with the issues trying to obtain a continuous distribution
    ## when dividing by an effectively discrete denominator
    splits <- seq(1, length(breaks) - 2)
    lmFits <- lapply(splits, function(x){
        vals <- seq(x, x + 2)
        fit <- lm.fit(x = cbind(1, breaks[vals]), y = freq[vals])
        fit$coefficients
    })

    ## Setup the values to return
    props <- seq(0, by = 1, length.out = bins) / (bins - 1)
    df <- tibble(
        GC_Content = props,
        Freq = vapply(props, function(x){
            interval <- findInterval(x, breaks[splits])
            max(lmFits[[interval]]["x2"]*x + lmFits[[interval]]["x1"], 0)
        }, numeric(1)))

    df

})