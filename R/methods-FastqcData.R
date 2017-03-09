#' Methods for an object of class FastqcData
#'
#' @description Methods and Accessors for an object of class FastqcData
#'
#' @details Descrbes all methods for an object of class \code{FastqcData}
#'
#' @param object An object of class \code{FastqcData}
#'
#' @include AllClasses.R
#'
#' @docType methods
#'
#' @importFrom scales comma
#'
#' @export
#' @rdname FastqcData-methods
#' @aliases path,FastqcData-method
setMethod("path", "FastqcData", function(object){object@path})

#' @export
#' @rdname FastqcData-methods
#' @aliases Basic_Statistics,FastqcData-method
setMethod("Basic_Statistics", "FastqcData", function(object){object@Basic_Statistics})

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_base_sequence_quality,FastqcData-method
setMethod("Per_base_sequence_quality", "FastqcData", function(object){object@Per_base_sequence_quality})

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_tile_sequence_quality,FastqcData-method
setMethod("Per_tile_sequence_quality", "FastqcData", function(object){object@Per_tile_sequence_quality})

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_sequence_quality_scores,FastqcData-method
setMethod("Per_sequence_quality_scores", "FastqcData", function(object){object@Per_sequence_quality_scores})

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_base_sequence_content,FastqcData-method
setMethod("Per_base_sequence_content", "FastqcData", function(object){object@Per_base_sequence_content})

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_sequence_GC_content,FastqcData-method
setMethod("Per_sequence_GC_content", "FastqcData", function(object){object@Per_sequence_GC_content})

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_base_N_content,FastqcData-method
setMethod("Per_base_N_content", "FastqcData", function(object){object@Per_base_N_content})

#' @export
#' @rdname FastqcData-methods
#' @aliases Sequence_Length_Distribution,FastqcData-method
setMethod("Sequence_Length_Distribution", "FastqcData", function(object){object@Sequence_Length_Distribution})

#' @export
#' @rdname FastqcData-methods
#' @aliases Sequence_Duplication_Levels,FastqcData-method
setMethod("Sequence_Duplication_Levels", "FastqcData", function(object){object@Sequence_Duplication_Levels})

#' @export
#' @rdname FastqcData-methods
#' @aliases Overrepresented_sequences,FastqcData-method
setMethod("Overrepresented_sequences", "FastqcData", function(object){object@Overrepresented_sequences})

#' @export
#' @rdname FastqcData-methods
#' @aliases Adapter_Content,FastqcData-method
setMethod("Adapter_Content", "FastqcData", function(object){object@Adapter_Content})

#' @export
#' @rdname FastqcData-methods
#' @aliases Kmer_Content,FastqcData-method
setMethod("Kmer_Content", "FastqcData", function(object){object@Kmer_Content})

#' @export
#' @rdname FastqcData-methods
#' @aliases Total_Deduplicated_Percentage,FastqcData-method
setMethod("Total_Deduplicated_Percentage", "FastqcData", function(object){object@Total_Deduplicated_Percentage})

#' @export
#' @rdname FastqcData-methods
#' @aliases Version,FastqcData-method
setMethod("Version", "FastqcData", function(object){object@Version})

#' @export
#' @rdname FastqcData-methods
#' @aliases Version,FastqcData-method
setMethod("getSummary", "FastqcData", function(object){object@Summary})

setMethod("show", "FastqcData",
          function(object){
            cat("FastqcData for", object@Basic_Statistics$Filename, "\n")
            cat("Source Fastq file contains", scales::comma(object@Basic_Statistics$Total_Sequences), "reads.\n")
            cat("Source FastQC file is located in", object@path)
          })
