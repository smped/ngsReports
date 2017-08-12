#' @title Return the File Names from an Object
#'
#' @description The File Names from an object
#'
#' @param object An object of class FastqcFile, FastqcFileList, FastqcData or FastqcDataList
#'
#' @return Returns the filenames, without the preceding directories, i.e. basename
#'
#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @name fileName
#' @aliases fileName,FastqcFile-method
#' @rdname fileName-methods
setMethod("fileName", "FastqcFile", function(object){basename(object@path)})

#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @name fileName
#' @aliases fileName,FastqcFileList-method
#' @rdname fileName-methods
setMethod("fileName", "FastqcFileList", function(object){vapply(object, fileName, character(1))})

#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @name fileName
#' @aliases fileName,FastqcData-method
#' @rdname fileName-methods
setMethod("fileName", "FastqcData", function(object){object@Summary$Filename[1]})

#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @name fileName
#' @aliases fileName,FastqcDataList-method
#' @rdname fileName-methods
setMethod("fileName", "FastqcDataList", function(object){vapply(object@.Data, fileName, character(1))})

