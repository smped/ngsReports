#' @title Return the File Names from an Object
#'
#' @description The File Names from an object
#'
#' @param object An object of class FastqcFile, FastqcFileList, FastqcData or FastqcDataList
#'
#' @return Returns the filenames, without the preceding directories, i.e. basename.
#' For a FastqcFile/FastqcFileList these will be the underlying FastQC files.
#' For a FastqcData/FastqcDataList, these will correspond to the fastq files 
#' about which the FastQC report was written.
#' 
#' @examples 
#' 
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#' 
#' # Form a FastqcFileList
#' ffl <- FastqcFileList(fileList)
#' fileName(ffl)
#' 
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' fileName(fdl)
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

