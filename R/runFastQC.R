#' Deprecated wrapper for the bash shell command fastqc.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' @details
#' `r lifecycle::badge("deprecated")`
#'
#' @param object Deprecated
#'
#' @export
#' @rdname runFastQC-methods
setGeneric("runFastQC", function(object){
    standardGeneric("runFastQC")
}
)
#' @aliases runFastQC,ANY-method
#' @rdname runFastQC-methods
#' @export
setMethod("runFastQC", "ANY", function(object){

    warning(
        "Deprecated. Please call `fastqc` directly using system2()."
    )

})
