#' Creates a Summarized Experiment from an accessibility matrix
#' 
#' Converts an accessibility matrix into a SummarizedExperiment file.
#' @param accessibility_matrix A matrix with first 3 columns containing coordinates. Remaining columns correspond to cells. Rows correspond to DNA sequences. Accessibility data should be binarized.
#' @param genome The genome the accessibility data was generated from (e.g. mm10, hg19).
#' @export
#' @examples 
#' matrix_to_summarized_experiment()

matrix_to_summarized_experiment<-function(accessibility_matrix,genome){
  coordinates<-accessibility_matrix[,1:3]
  data<-as.matrix(accessibility_matrix[,-c(1:3)])
  assay<-SummarizedExperiment(data)
  rowRanges(assay)<-GRanges(coordinates)
  assayNames(assay)<-"counts"
  assay<-add_gcbias(assay,genome = genome)
}