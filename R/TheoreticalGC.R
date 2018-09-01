#' @title The TheoreticalGC Object Class
#'
#' @description Contains Theoretical GC content for a selection of species
#'
#' @details Estimates are able to be retained for genomic and transcriptomic sequences
#' Values are stored as frequencies
#'
#' @return An object of class TheoreticalGC
#'
#' @include validationFunctions.R
#'
#' @slot Genome A \code{data.frame} containing theoretical GC content for genomic sequences
#' @slot Transcriptome A \code{data.frame} containing theoretical GC content for transcriptomic sequences
#' @slot mData A \code{data.frame} containing metadata about all species in the object
setClass("TheoreticalGC", slots=c(Genome = "data.frame",
                                  Transcriptome = "data.frame",
                                  mData = "data.frame"))
setValidity("TheoreticalGC", isValidTheoreticalGC)

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

#' @title List Available Genomes
#'
#' @description List available genomes in a TheoreticalGC object
#'
#' @param object An object of class Theoretical GC
#'
#' @return A \code{tibble} object
#' 
#' @examples 
#' genomes(gcTheoretical)
#'
#' @export
#' @name genomes
#' @rdname genomes
setGeneric("genomes", function(object){standardGeneric("genomes")})


#' @importFrom methods slot
#' @export
#' @rdname genomes
#' @aliases genomes,TheoreticalGC-method
setMethod("genomes", "TheoreticalGC", function(object){
    gn <- object@mData$Genome
    dplyr::select(object@mData[gn,], -tidyselect::ends_with("ome"))
})

#' @title List Available Transcriptomes
#'
#' @description List available transcriptomes in a TheoreticalGC object
#'
#' @param object An object of class Theoretical GC
#'
#' @return A \code{tibble} object
#'
#' @examples 
#' transcriptomes(gcTheoretical)
#'
#' @export
#' @name transcriptomes
#' @rdname transcriptomes
setGeneric("transcriptomes", function(object){standardGeneric("transcriptomes")})

#' @importFrom methods slot
#' @export
#' @rdname transcriptomes
#' @aliases transcriptomes,TheoreticalGC-method
setMethod("transcriptomes", "TheoreticalGC", function(object){
    tr <- object@mData$Transcriptome
    dplyr::select(object@mData[tr,], -tidyselect::ends_with("ome"))
})

#' @title Get Theoretical GC content
#'
#' @description Get the GC content data from a TheoreticalGC object
#'
#' @param object An object of class Theoretical GC
#' @param name The Name of the species in 'Gspecies' format, e.g. Hsapiens
#' @param type The type of GC content. Can only be either "Genome" or "Transcriptome"
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
setMethod("getGC", "TheoreticalGC",
          function(object, name, type){
              
              type <- stringr::str_to_title(type)
              type <- match.arg(type[[1]], c("Genome", "Transcriptome"))
              
              if (type == "Genome"){
                  col <- tryCatch(match.arg(name, colnames(object@Genome)))
                  df <- object@Genome[c("GC_Content", col)]
              }
              else{
                  col <- tryCatch(match.arg(name, colnames(object@Transcriptome)))
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
