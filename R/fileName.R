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
setMethod("fileName", "FastqcFile", function(object){basename(object@path)})

#' @export
setMethod("fileName", "FastqcFileList", function(object){vapply(object, fileName, character(1))})

#' @export
setMethod("fileName", "FastqcData", function(object){object@Summary$Filename[1]})

#' @export
setMethod("fileName", "FastqcDataList", function(object){vapply(object@.Data, fileName, character(1))})

