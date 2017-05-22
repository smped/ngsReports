# Define object classes
setClass("FastqcFile", slots = c(path = "character"))
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
setClass("PwfCols", slots = c(PASS = "character",
                              WARN = "character",
                              FAIL = "character",
                              MAX = "character"))

# Set the validation functions for any object classes
#' @include validationFunctions.R
setValidity("FastqcFile", isValidFastqcFile)
setValidity("FastqcFileList", isValidFastqcFileList)
# setValidity("FastqcData", isValidFastqcData) # Not written or defined yet
setValidity("FastqcDataList", isValidFastqcDataList) # Not written or defined yet
setValidity("PwfCols", isValidPwf)

# These are never set
setMethod("names<-", "FastqcFile", function(x){
  warning("The names attribute cannot be set on a FastqcFile object")
  x
  })
setMethod("names<-", "FastqcFileList", function(x){
  warning("The names attribute cannot be set on a FastqcFileList object")
  x
})
setMethod("names<-", "FastqcData", function(x){
  warning("The names attribute cannot be set on a FastqcData object")
  x
})
setMethod("names<-", "FastqcDataList", function(x){
  warning("The names attribute cannot be set on a FastqcDataList object")
  x
})
setMethod("names<-", "PwfCols", function(x){
  warning("The names attribute cannot be set on a PwfCols object")
  x
})
