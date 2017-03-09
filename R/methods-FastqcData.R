#' Methods for an object of class FastqcData
#'
#' @description Methods and Accessors for an object of class FastqcData
#'
#' @details Describes all methods for an object of class \code{FastqcData}
#'
#' @param object An object of class \code{FastqcData}
#'
#' @include AllClasses.R
#'
#' @docType methods
#'
#' @importFrom scales comma
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr everything
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
setMethod("Per_base_sequence_quality", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_base_sequence_quality,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
            })

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_tile_sequence_quality,FastqcData-method
setMethod("Per_tile_sequence_quality", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_tile_sequence_quality,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
            })

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_sequence_quality_scores,FastqcData-method
setMethod("Per_sequence_quality_scores", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_sequence_quality_scores,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
            })

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_base_sequence_content,FastqcData-method
setMethod("Per_base_sequence_content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_base_sequence_content,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_sequence_GC_content,FastqcData-method
setMethod("Per_sequence_GC_content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_sequence_GC_content,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Per_base_N_content,FastqcData-method
setMethod("Per_base_N_content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_base_N_content,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Sequence_Length_Distribution,FastqcData-method
setMethod("Sequence_Length_Distribution", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Sequence_Length_Distribution,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Sequence_Duplication_Levels,FastqcData-method
setMethod("Sequence_Duplication_Levels", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Sequence_Duplication_Levels,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Overrepresented_sequences,FastqcData-method
setMethod("Overrepresented_sequences", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Overrepresented_sequences,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Adapter_Content,FastqcData-method
setMethod("Adapter_Content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Adapter_Content,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Kmer_Content,FastqcData-method
setMethod("Kmer_Content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Kmer_Content,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Total_Deduplicated_Percentage,FastqcData-method
setMethod("Total_Deduplicated_Percentage", "FastqcData",
          function(object){
            dplyr::data_frame(Filename = fileNames(object),
                              Total = object@Total_Deduplicated_Percentage)
          })

#' @export
#' @rdname FastqcData-methods
#' @aliases Version,FastqcData-method
setMethod("Version", "FastqcData", function(object){object@Version})

#' @export
#' @rdname FastqcData-methods
#' @aliases getSummary,FastqcData-method
setMethod("getSummary", "FastqcData", function(object){object@Summary})

#' @export
#' @rdname FastqcData-methods
#' @aliases names,FastqcData-method
setMethod("names", "FastqcData", function(x){basename(x@path)})

#' @export
#' @rdname FastqcData-methods
#' @aliases fileName,FastqcData-method
setMethod("fileNames", "FastqcData", function(object){object@Summary$FastqFile[1]})

setMethod("show", "FastqcData",
          function(object){
            cat("FastqcData for", object@Basic_Statistics$Filename, "\n")
            cat("Source Fastq file contains", scales::comma(object@Basic_Statistics$Total_Sequences), "reads.\n")
            cat("Source FastQC file is located in", object@path)
          })
