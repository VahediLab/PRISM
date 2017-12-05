#' Bedgraph files to SummarizedExperiment
#' 
#' Creates a SummarizedExperiment from a list of bedgraph files for later analyses.
#' @param directory_bedgraph Directory containing a list of bedgraph files. First 3 columns correspond to coordinates. User can specify which column contains accessibility values (automatically column 4).
#' @param column_accessibility The column containing the number of reads. Default is column 4.
#' @param genome The genome the accessibility data was generated from (e.g. mm10, hg19).
#' @export
#' @examples 
#' bedgraphs_to_summarized_experiment()

bedgraphs_to_summarized_experiment<-function(directory_bedgraph,genome, column_accessibility = 4){
  accessibility_matrix<-bedgraphs_to_accessibility_matrix(directory_bedgraph,column_accessibility)
  matrix_to_summarized_experiment(accessibility_matrix,genome)
}
