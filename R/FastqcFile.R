#' Create a New FastqcFile Object
#'
#' Create a FastqcFile object
#'
#' @param x Character vector (1) specifying a valid path to a file/directory as output by FastQC
#'
#' @docType methods
setGeneric("FastqcFile",function(x){standardGeneric("FastqcFile")})

#' @importFrom methods new
#' @export
setMethod("FastqcFile", "character",
          function(x){
            new("FastqcFile", path = x)
          })


