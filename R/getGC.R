#' @title Get Theoretical GC content
#'
#' @description Get the GC content data from a TheoreticalGC object
#'
#' @param object An object of class Theoretical GC
#' @param name The Name of the species in 'Gspecies' format, e.g. Hsapiens
#' @param type The type of GC content. Can only be either "Genome" or
#' "Transcriptome"
#'
#' @return A `tibble` object
#'
#' @examples
#' getGC(gcTheoretical, name = "Hsapiens", type = "Genome")
#'
#' @docType methods
#'
#' @export
#' @name getGC
#' @rdname getGC
setGeneric("getGC", function(object, name, type){standardGeneric("getGC")})
#' @export
#' @rdname getGC
setMethod("getGC", "ANY", function(object, type){
    .errNotImp(object)
})
#' @importFrom methods slot
#' @rdname getGC
#' @export
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
