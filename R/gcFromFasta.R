#' @title Calculate GC content from a fasta file
#' 
#' @description Calculate GC content from a set of DNA sequences provided by a fasta file
#' 
#' @details Used to obtain theoretical GC content for a given genome/transcriptome
#' 
#' @param Fastafile A fasta file containing a genome or transcriptome of interest
#' @param n the number of fragments to simulate
#' @param bp the fragment length used in estimation
#' 
#' @return A data.frame with the GC percentage and the exected frequency
#' 
#' @author Hien To <totuhien@gmail.com>
#' 
#' @keywords internal
gcFromFasta <- function(Fastafile, n=1e+6, bp=100){
    ref <- readDNAStringSet(filepath = Fastafile)
    file.gc = "gc.txt"
    fastqcTheoreticalGC::generateDistn(ref, file = file.gc, n = n, bp = bp)
    gc <- data.frame(read.table(file.gc, header = FALSE, stringsAsFactors = FALSE))
    file.remove(file.gc)
    colnames(gc) <- c("GC_Content","Freq")
    gc$Freq <- gc$Freq/sum(gc$Freq)
    gc
}