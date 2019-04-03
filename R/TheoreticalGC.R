#' @title The TheoreticalGC Object Class
#'
#' @description Contains Theoretical GC content for a selection of species
#'
#' @details Estimates are able to be retained for genomic and transcriptomic
#' sequences. Values are stored as frequencies.
#'
#' @return An object of class TheoreticalGC
#'
#' @include validationFunctions.R
#'
#' @slot Genome A \code{data.frame} containing theoretical GC content for
#' genomic sequences
#' @slot Transcriptome A \code{data.frame} containing theoretical GC content
#' for transcriptomic sequences
#' @slot mData A \code{data.frame} containing metadata about all species in the
#' object
#'
#' @examples
#'
#' ## How to form an object using your own fasta file
#' faDir <- system.file("extdata", package = "ngsReports")
#' faFile <- list.files(faDir, pattern = "fasta", full.names = TRUE)
#' gen_df <- getGcDistn(faFile, n = 200)
#' gen_df <- dplyr::rename(gen_df, Athaliana = Freq)
#' mData_df <-
#'     data.frame(Name = "Athaliana", Genome = TRUE, Transcriptome = FALSE)
#' tr_df <- data.frame()
#' myGC <- new(
#'    "TheoreticalGC", Genome = gen_df, Transcriptome = tr_df, mData = mData_df)
#'
setClass("TheoreticalGC",  slots = c(
    Genome = "data.frame",
    Transcriptome = "data.frame",
    mData = "data.frame")
)
setValidity("TheoreticalGC", .isValidTheoreticalGC)

#' @title Extract Metadata for TheoreticalGC objects
#'
#' @description Extract Metadata for TheoreticalGC objects
#'
#' @param object An object of class Theoretical GC
#'
#' @return A \code{tibble} object
#'
#' @examples
#' mData(gcTheoretical)
#'
#' @export
#' @name mData
#' @rdname mData
setGeneric("mData", function(object){standardGeneric("mData")})

#' @importFrom methods slot
#' @export
#' @rdname mData
#' @aliases mData,TheoreticalGC-method
setMethod("mData", "TheoreticalGC", function(object){object@mData})

#' @title List Genomes or Transcriptomes with Theoretical GC Content
#'
#' @description
#' List available genomes or transcriptomes in a TheoreticalGC object
#'
#' @details
#' An object of class TheoreticalGC can hold the theoretical GC content for one
#' or more species, for either the genome or transriptome.
#' This function checks which species are available in the given object, for
#' either the genome or transcriptome, as supplied to the parameter \code{type}.
#'
#' @param object An object of class TheoreticalGC
#' @param type character indicating either Genome or Transcriptome
#'
#' @return A \code{tibble} object
#'
#' @examples
#' gcAvail(gcTheoretical, "Genome")
#'
#' @export
setGeneric("gcAvail", function(object, type){standardGeneric("gcAvail")})


#' @importFrom methods slot
#' @export
#' @rdname gcAvail
#' @aliases gcAvail,TheoreticalGC-method
setMethod(
    "gcAvail",
    "TheoreticalGC",
    function(object, type){
        type <- match.arg(type, c("Genome", "Transcriptome"))
        ln <- slot(object, "mData")[[type]]
        df <- slot(object, "mData")[ln,]
        dplyr::select(df, -tidyselect::ends_with("ome"))
    })

#' @title Get Theoretical GC content
#'
#' @description Get the GC content data from a TheoreticalGC object
#'
#' @param object An object of class Theoretical GC
#' @param name The Name of the species in 'Gspecies' format, e.g. Hsapiens
#' @param type The type of GC content. Can only be either "Genome" or
#' "Transcriptome"
#'
#' @return A \code{tibble} object
#'
#' @examples
#' getGC(gcTheoretical, name = "Hsapiens", type = "Genome")
#'
#' @export
#' @name getGC
#' @rdname getGC
setGeneric("getGC", function(object, name, type){standardGeneric("getGC")})

#' @importFrom methods slot
#' @export
#' @rdname getGC
#' @aliases getGC,TheoreticalGC-method
setMethod("getGC", "TheoreticalGC", function(object, name, type){

    type <- stringr::str_to_title(type)
    type <- match.arg(type[[1]], c("Genome", "Transcriptome"))

    if (type == "Genome") {
        col <- match.arg(name, colnames(object@Genome))
        df <- object@Genome[c("GC_Content", col)]
    }
    if (type == "Transcriptome") {
        col <- match.arg(name, colnames(object@Transcriptome))
        df <- object@Transcriptome[c("GC_Content", col)]
    }

    df

})


setMethod("show", "TheoreticalGC", function(object){
    meta <- mData(object)
    nGenomes <- sum(meta$Genome)
    nTranscriptomes <- sum(meta$Transcriptome)
    cat("TheoreticalGC Object for:\n")
    cat(nGenomes, "Genomes &", nTranscriptomes, "Transcriptomes\n")
})
