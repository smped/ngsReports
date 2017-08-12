# Define object classes
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


# Set the validation functions for any object classes
#' @include validationFunctions.R
# setValidity("FastqcData", isValidFastqcData) # Not written or defined yet
