#' Calculate unnormalized cell-to-cell variation in chromatin accessibility
#' 
#' Calculates unnormalized cell-to-cell variation in chromatin accessibility using angular cosine distance inputted into principal coordinate analysis. Then each cell's distance from the centroid comprises the accessibility variation. Technical biases associated with GC content and mean accessibility are controlled for.
#' @param scATACseq_experiment A SummarizedExperiment containing scATACseq data.
#' @param ChIPseq_experiment A matrix containing ChIPseq data for matching DNA sequences.
#' @examples 
#' raw_variation()

raw_variation<-function(scATACseq_experiment,ChIPseq_experiment){
  treatment.rows<-which(ChIPseq_experiment==1)
  control.rows<-which(ChIPseq_experiment==0)
  treatmentMatrix<-assay(scATACseq_experiment[treatment.rows,])
  controlMatrix<-assay(scATACseq_experiment[control.rows,])
  answers<-get_background_peaks(scATACseq_experiment,ChIPseq_experiment)
  Background<-answers[[1]]
  eligibleMatrix<-answers[[2]]
  rm<-which(rowSums(Background)==0)
  if(isEmpty(rm)){
    Back_Variability<-mean(unlist(lapply(1:30,function(j){
      mean(centroid_distance(eligibleMatrix[Background[,j],]))})))
    Front_Variability<-mean(centroid_distance(treatmentMatrix))
  } else {
    Back_Variability<-mean(unlist(lapply(1:30,function(j){
      mean(centroid_distance(eligibleMatrix[Background[,j],][-rm,]))})))
    Front_Variability<-mean(centroid_distance(treatmentMatrix[-rm,]))
  }
  return(Front_Variability/Back_Variability)
}