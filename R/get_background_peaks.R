#' Background peak selection
#' 
#' Selects background peaks to control for technical biases associated with GC content and mean accessibility. Background peaks are drawn from non-TF bound peaks and matched peak-for-peak for GC content, and matrix-for-matrix for mean accessibility.
#' @param scATACseq_experiment A SummarizedExperiment with scATACseq data.
#' @param ChIPseq_experiment A ChIPseq matrix.
#' @export
#' @examples 
#' get_background_peaks()

get_background_peaks<-function(scATACseq_experiment,ChIPseq_experiment){
  treatment.rows<-which(ChIPseq_experiment==1)
  control.rows<-which(ChIPseq_experiment==0)
  treatmentMatrix<-assay(scATACseq_experiment[treatment.rows,])
  controlMatrix<-assay(scATACseq_experiment[control.rows,])
  meanAccessibility<-mean(treatmentMatrix)
  mean.matched.rows<-which(rowMeans(controlMatrix)<(meanAccessibility+0.01)&rowMeans(controlMatrix)>(meanAccessibility-0.01))
  eligibleMatrix<-controlMatrix[mean.matched.rows,]
  gcData<-rowData(scATACseq_experiment[mean.matched.rows,])$bias
  gcQuantiles<-quantile(gcData,seq(0,1,0.02))
  gcControl<-lapply(2:51,function(i){which(gcData<gcQuantiles[[i]]&gcData>=gcQuantiles[[i-1]])})
  gcTreatment<-findInterval(rowData(scATACseq_experiment[treatment.rows,])$bias,gcQuantiles)
  back<-lapply(1:length(treatment.rows),function(i){
    if (length(unlist(gcControl[gcTreatment[[i]]]))==0) {replicate(30,0)}
    else {sample(unlist(gcControl[gcTreatment[[i]]]),30,replace=TRUE)}
  })
  Background<-do.call(rbind,back)
  return(list(Background,eligibleMatrix))
}