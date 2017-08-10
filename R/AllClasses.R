# Define object classes
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

setClass("TheoreticalGC", slots=c(Alyrata = "data.frame",
                                  Amellifera = "data.frame",
                                  Athaliana = "data.frame",
                                  Btaurus = "data.frame",
                                  Celegans = "data.frame",
                                  Cfamiliaris = "data.frame",
                                  Dmelanogaster = "data.frame",
                                  Drerio = "data.frame",
                                  Ecoli = "data.frame",
                                  Gaculeatus = "data.frame",
                                  Ggallus = "data.frame",
                                  Hsapiens = "data.frame",
                                  Mfascicularis = "data.frame",
                                  Mfuro = "data.frame",
                                  Mmulatta = "data.frame",
                                  Mmusculus = "data.frame",
                                  Osativa = "data.frame",
                                  Ptroglodytes = "data.frame",
                                  Rnorvegicus = "data.frame",
                                  Scerevisiae = "data.frame",
                                  Sscrofa = "data.frame",
                                  Tgondii = "data.frame",
                                  Tguttata = "data.frame",
                                  Vvinifera = "data.frame",
                                  mData = "data.frame"))

# Set the validation functions for any object classes
#' @include validationFunctions.R
setValidity("FastqcFileList", isValidFastqcFileList)
# setValidity("FastqcData", isValidFastqcData) # Not written or defined yet
setValidity("FastqcDataList", isValidFastqcDataList) # Not written or defined yet
setValidity("TheoreticalGC", isValidTheoreticalGC)

# These are never set
setMethod("names<-", "FastqcFile", function(x, value){
  warning("The names attribute cannot be set on a FastqcFile object")
  x
})
setMethod("names<-", "FastqcFileList", function(x, value){
  warning("The names attribute cannot be set on a FastqcFileList object")
  x
})
setMethod("names<-", "FastqcData", function(x, value){
  warning("The names attribute cannot be set on a FastqcData object")
  x
})
setMethod("names<-", "FastqcDataList", function(x, value){
  warning("The names attribute cannot be set on a FastqcDataList object")
  x
})
