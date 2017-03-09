#' Methods for an object of class FastqcDataList
#'
#' @description Methods and Accessors for an object of class FastqcDataList
#'
#' @details Describes all methods for an object of class \code{FastqcDataList}
#'
#' @param object An object of class \code{FastqcDataList}
#'
#' @return Returns a \code{data_frame} for all the methods except \code{path}, \code{names},
#' \code{fileNames}
#'
#' @include AllClasses.R
#'
#' @docType methods
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#' @rdname FastqcDataList-methods
#' @aliases path,FastqcDataList-methods
setMethod("path", "FastqcDataList", function(object){vapply(object@.Data, path, character(1))})

#' @rdname FastqcDataList-methods
#' @aliases names,FastqcDataList-methods
setMethod("names", "FastqcDataList", function(x){vapply(x@.Data, names, character(1))})

#' @export
#' @rdname FastqcDataList-methods
#' @aliases fileName,FastqcDataList-method
setMethod("fileNames", "FastqcDataList", function(object){vapply(object@.Data, fileNames, character(1))})

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Basic_Statistics,FastqcDataList-method
setMethod("Basic_Statistics", "FastqcDataList",
          function(object){
            x <- lapply(object@.Data, Basic_Statistics)
            dplyr::bind_rows(x)
            })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Per_base_sequence_quality,FastqcDataList-method
setMethod("Per_base_sequence_quality", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_base_sequence_quality)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Per_tile_sequence_quality,FastqcDataList-method
setMethod("Per_tile_sequence_quality", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_tile_sequence_quality)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Per_sequence_quality_scores,FastqcDataList-method
setMethod("Per_sequence_quality_scores", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_sequence_quality_scores)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Per_base_sequence_content,FastqcDataList-method
setMethod("Per_base_sequence_content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_base_sequence_content)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Per_sequence_GC_content,FastqcDataList-method
setMethod("Per_sequence_GC_content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_sequence_GC_content)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Per_base_N_content,FastqcDataList-method
setMethod("Per_base_N_content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_base_N_content)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Sequence_Length_Distribution,FastqcDataList-method
setMethod("Sequence_Length_Distribution", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Sequence_Length_Distribution)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Sequence_Duplication_Levels,FastqcDataList-method
setMethod("Sequence_Duplication_Levels", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Sequence_Duplication_Levels)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Overrepresented_sequences,FastqcDataList-method
setMethod("Overrepresented_sequences", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Overrepresented_sequences)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Adapter_Content,FastqcDataList-method
setMethod("Adapter_Content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Adapter_Content)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Kmer_Content,FastqcDataList-method
setMethod("Kmer_Content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Kmer_Content)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Total_Deduplicated_Percentage,FastqcDataList-method
setMethod("Total_Deduplicated_Percentage", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Total_Deduplicated_Percentage)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases getSummary,FastqcDataList-method
setMethod("getSummary", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, getSummary)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname FastqcDataList-methods
#' @aliases Version,FastqcDataList-method
setMethod("Version", "FastqcDataList",
          function(object){
            data_frame(Filename = fileNames(object),
                       Version = vapply(object@.Data, Version, character(1)))
          })

setMethod("show", "FastqcDataList",
          function(object){
            l <- length(object)
            cat("FastqcDataList for", l, "file(s).\n")
          })

