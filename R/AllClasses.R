# Define object classes
setClass("FastqcFile", slots = c(path = "character", compressed = "logical"))
setClass("FastqcFileList", contains="list")
setClass("FastqcData", slots = c(Summary = "data.frame",
                                 Basic_Statistics = "data.frame",
                                 Per_base_sequence_quality = "data.frame",
                                 Per_tile_sequence_quality = "data.frame",
                                 Per_sequence_quality_scores = "data.frame",
                                 Per_base_sequence_content = "data.frame",
                                 Per_sequence_GC_content = "data.frame",
                                 Per_base_N_content = "data.frame",
                                 Sequence_Length_Distribution = "data.frame",
                                 Sequence_Duplication_Levels = "data.frame",
                                 Overrepresented_sequences = "data.frame",
                                 Adapter_Content = "data.frame",
                                 Kmer_Content = "data.frame",
                                 Total_Deduplicated_Percentage = "numeric",
                                 Version = "character",
                                 path = "character"))
setClass("FastqcDataList", contains="list")

# Set the vailation functions for any object classes
#' @include validationFunctions.R
setValidity("FastqcFile", validFastqcFile)
setValidity("FastqcFileList", validFastqcFileList)
# setValidity("FastqcData", validFastqcData) # Not written or defined yet
# setValidity("FastqcDataList", validFastqcData) # Not written or defined yet

# Set the Generics
setGeneric("FastqcFile",function(filePath){standardGeneric("FastqcFile")})
setGeneric("FastqcFileList", function(path){standardGeneric("FastqcFileList")})
setGeneric("isCompressed", function(object){standardGeneric("isCompressed")})
setGeneric("path", function(object){standardGeneric("path")})
setGeneric("getSummary", function(object){standardGeneric("getSummary")})
setGeneric("getFastqcData", function(object){standardGeneric("getFastqcData")})
setGeneric("Basic_Statistics", function(object){standardGeneric("Basic_Statistics")})
setGeneric("Per_base_sequence_quality", function(object){standardGeneric("Per_base_sequence_quality")})
setGeneric("Per_tile_sequence_quality", function(object){standardGeneric("Per_tile_sequence_quality")})
setGeneric("Per_sequence_quality_scores", function(object){standardGeneric("Per_sequence_quality_scores")})
setGeneric("Per_base_sequence_content", function(object){standardGeneric("Per_base_sequence_content")})
setGeneric("Per_sequence_GC_content", function(object){standardGeneric("Per_sequence_GC_content")})
setGeneric("Per_base_N_content", function(object){standardGeneric("Per_base_N_content")})
setGeneric("Sequence_Length_Distribution", function(object){standardGeneric("Sequence_Length_Distribution")})
setGeneric("Sequence_Duplication_Levels", function(object){standardGeneric("Sequence_Duplication_Levels")})
setGeneric("Overrepresented_sequences", function(object){standardGeneric("Overrepresented_sequences")})
setGeneric("Adapter_Content", function(object){standardGeneric("Adapter_Content")})
setGeneric("Kmer_Content", function(object){standardGeneric("Kmer_Content")})
setGeneric("Total_Deduplicated_Percentage", function(object){standardGeneric("Total_Deduplicated_Percentage")})
setGeneric("Version", function(object){standardGeneric("Version")})
setGeneric("fileNames",function(object){standardGeneric("fileNames")})

