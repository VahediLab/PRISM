#' Calculate GC content of every DNA sequence
#' 
#' Calculates GC content of every DNA sequence from a SummarizedExperiment file. 
#' @param object A SummarizedExperiment file
#' @param genome The genome the accessibility data was generated from (e.g. mm10, hg19).
#' @export
#' @examples 
#' add_gcbias()

add_gcbias<-function(object, 
                   genome = GenomeInfoDb::genome(object)) {
            peaks <- SummarizedExperiment::rowRanges(object)
            seqs <- BSgenome::getSeq(genome, peaks)
            nucfreqs <- BSgenome::letterFrequency(seqs, c("G", "C", "A", "T"))
            gc <- apply(nucfreqs, 1, function(x) rowSums(x[,1:2])/rowSums(x))
            SummarizedExperiment::rowRanges(object)$bias <- gc
            return(object)
}