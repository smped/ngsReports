#' @title Return the FileNames from an Object
#'
#' @description The File Names from an object
#'
#' @param object An object of class FastqcFile, FastqcFileList, FastqcData or FastqcDataList
#'
#' @return Returns the filenames, without the preceding directories, i.e. basename
#'
setGeneric("fileNames",function(object){standardGeneric("fileNames")})

#' @export
setMethod("fileNames", "FastqcFile", function(object){basename(object@path)})

#' @export
setMethod("fileNames", "FastqcFileList", function(object){vapply(object, fileNames, character(1))})

#' @export
setMethod("fileNames", "FastqcData", function(object){object@Summary$Filename[1]})

#' @export
setMethod("fileNames", "FastqcDataList", function(object){vapply(object@.Data, fileNames, character(1))})

