#' @title Theoretical GC content
#'
#' @description This object contains the theoretical GC content for each
#' provided species, for both the genome and transcriptome, where available.
#'
#' @details The object is defined with the S4 class \code{TheoreticalGC}.
#' Species for which information is available can be found using
#' the commands \code{genomes(gcTheoretical)}
#' or \code{transcriptomes(gcTheoretical)}.
#'
#' Metadata is accessible using \code{mData(gcTheoretical)}.
#'
#' All GC content was calculated using code from
#' https://github.com/mikelove/fastqcTheoreticalGC using BSgenome packages.
#' This provides a default set of GC content data for common organisms
#' generated using 100bp reads/fragments and 1e6 reads.
#'
#' @examples
#' ## Check which genomes are included
#' genomes(gcTheoretical)
#'
#' ## Check which transcriptomes are included
#' transcriptomes(gcTheoretical)
#'
#' ## Metadata is also included
#' mData(gcTheoretical)
#'
"gcTheoretical"
