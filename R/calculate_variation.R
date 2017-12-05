#' Calculate cell-to-cell variation in chromatin accessibility
#' 
#' Calculates unnormalized cell-to-cell variation in chromatin accessibility using angular cosine distance inputted into principal coordinate analysis. Then each cell's distance from the centroid comprises the accessibility variation. Technical biases associated with GC content and mean accessibility are controlled for.
#' @param scATACseq_experiment A SummarizedExperiment containing scATACseq data.
#' @param ChIPseq_experiment A matrix containing ChIPseq data for matching DNA sequences.
#' @export
#' @examples 
#' calculate_variation()

calculate_variation<-function(scATACseq_experiment,ChIPseq_experiment){
  unnormalized_variation<-raw_variation(scATACseq_experiment,ChIPseq_experiment)
  treatment.rows<-which(ChIPseq_experiment==1)
  control.rows<-which(ChIPseq_experiment==0)
  lit<-mean(ChIPseq_experiment)
  fakes<-unlist(lapply(1:30,function(j){
    fake<-as.matrix(ifelse(runif(nrow(B1.1))<lit,1,0))
    raw_variation(B1.1,fake)
  }))
  normalizing_factor<-mean(fakes)
  normalized_variation<-unnormalized_variation/normalizing_factor
  Zscore<-(unnormalized_variation-normalizing_factor)/sd(fakes)
  pvalue<-2*pnorm(-abs(Zscore))
  return(c(normalized_variation, pvalue))
}