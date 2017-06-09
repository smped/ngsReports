#' Methods for an object of class TheoreticalGC
#'
#' @description Methods and Accessors for an object of class TheoreticalGC
#'
#' @details Describes all methods for an object of class \code{TheoreticalGC}.
#'
#' \code{getGC} returns the theoretical GC content for a given species as a data frame.
#'
#' \code{mData} returns the metadata data frame for all species.
#'
#' \code{genomes} and \code{transcriptomes} return vectors listing the species for which the
#' theoretical GC content is available.
#'
#' \code{allSpecies} returns a vector of the species for which the theoretical GC content is available
#'
#' @param object An object of class \code{TheoreticalGC}
#' @param species character. The species to extract data for
#' @param type character. Can be \code{"Genomic"} or \code{"Transcriptomic"}
#'
#' @return A \code{data_frame}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#'
#' @export
setMethod("getGC", "TheoreticalGC",
          function(object, species, type){

            species <- tryCatch(match.arg(species, slotNames(object)))
            df <- slot(object, species)

            type <- match.arg(type[1], colnames(slot(object, species)))
            if(!type %in% colnames(df)) stop(paste("No", type, "information available for", species))

            dplyr::select(df, GC_Content, one_of(type))

          })

#' @export
#' @rdname mData,TheoreticalGC-method
setMethod("mData", "TheoreticalGC", function(object){object@mData})

#' @export
#' @rdname allSpecies,TheoreticalGC-method
setMethod("allSpecies", "TheoreticalGC", function(object){object@mData$Species})

#' @export
#' @rdname genomes,TheoreticalGC-method
setMethod("genomes", "TheoreticalGC", function(object){
  dplyr::filter(object@mData, Genome)$Species
  })

#' @export
#' @rdname transcriptomes,TheoreticalGC-method
setMethod("transcriptomes", "TheoreticalGC", function(object){
  dplyr::filter(object@mData, Transcriptome)$Species
})


# The show method doesn't need exporting
setMethod("show", "TheoreticalGC",
          function(object){
            n <- nrow(object@mData)
            cat("Theoretical GC content for", n, "species")
          })
