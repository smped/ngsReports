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

# Set the validation functions for any object classes
#' @include validationFunctions.R
setValidity("FastqcFile", validFastqcFile)
setValidity("FastqcFileList", validFastqcFileList)
# setValidity("FastqcData", validFastqcData) # Not written or defined yet
# setValidity("FastqcDataList", validFastqcData) # Not written or defined yet

# Set the show methods that don't need exporting
setMethod("show", "FastqcDataList",
          function(object){
            l <- length(object)
            cat("FastqcDataList for", l, "file(s).\n")
          })
setMethod("show", "FastqcData",
          function(object){
            cat("FastqcData for", object@Basic_Statistics$Filename, "\n")
            cat("Source Fastq file contains", scales::comma(object@Basic_Statistics$Total_Sequences), "reads.\n")
            cat("Source FastQC file is located in", object@path)
          })
setMethod("show", "FastqcFileList",
          function(object){
            l <- length(object)
            cmp <- sum(isCompressed(object))
            cat("FastqcFileList of", l, "file(s).\n")
            cat("Located in:\n", paste(unique(dirname(path(object))), collapse = "\n"))
          })
setMethod("show", "FastqcFile",
          function(object){
            cat(fileNames(object), "\n")
            cat("Located in", dirname(path(object)), "\n")
          })

# These are never set
setMethod("names", "FastqcFile", NULL)
setMethod("names", "FastqcFileList", NULL)
setMethod("names<-", "FastqcFile", function(x){
  warning("The names attribute cannot be set on a FastqcFile object")
  x
  })
setMethod("names<-", "FastqcFileList", function(x){
  warning("The names attribute cannot be set on a FastqcFile object")
  x
})
