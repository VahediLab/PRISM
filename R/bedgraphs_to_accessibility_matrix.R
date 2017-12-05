#' Bedgraph files to accessibility matrix
#' 
#' Creates an accessibility matrix from a list of bedgraph files. First 3 columns are coordinates, remaining columns correspond to cells, and every row is a DNA sequence.
#' @param directory_bedgraph Directory containing a list of bedgraph files. First 3 columns correspond to coordinates. User can specify which column contains accessibility values (automatically column 4).
#' @param column_accessibility The column containing the number of reads. Default is column 4.
#' @export
#' @examples 
#' bedgraphs_to_accessibility_matrix()

bedgraphs_to_accessibility_matrix<-function(directory_bedgraph,column_accessibility = 4){
  setwd(directory_bedgraph)
  fileList <- list.files(path=directory_bedgraph, pattern=".bedgraph")
  tableList <- sapply(fileList, read.table)
  accessibilityList <- lapply(1:ncol(tableList),function(i){cbind(paste(tableList[,i][[1]],tableList[,i][[2]],tableList[,i][[3]]),tableList[,i][[column_accessibility]])})
  preCoordinates<-unique(unlist(lapply(1:length(accessibilityList),function(i){accessibilityList[[i]][,1]})))
  preCoordinates<-preCoordinates[order(preCoordinates)]
  processedList<-lapply(1:length(accessibilityList),function(i){
    accessibilityMatrix<-matrix(nrow=length(preCoordinates),ncol=2)
    accessibilityMatrix[1:nrow(unique(accessibilityList[[i]])),]<-unique(accessibilityList[[i]])
    return (accessibilityMatrix[order(accessibilityMatrix[,1]),])
  })
  processedMatrix<-do.call(cbind, processedList)
  pre_scATACseq<-processedMatrix[,c(1,seq(2,ncol(processedMatrix),2))]
  accessibilityData<-pre_scATACseq[,-1]
  mode(accessibilityData)<-"numeric"
  rm.accessibilityData<-which(rowSums(accessibilityData)==0)
  if(!isEmpty(rm.accessibilityData)) {
    accessibilityData2<-accessibilityData[-rm.accessibilityData,]
    binaryAccessibility<-ifelse(accessibilityData2>0,1,0)
    coordinates<-read.table(text = pre_scATACseq[,1], sep = " ", colClasses = "character")
    coordinates2<-coordinates[-rm.accessibilityData,]
    colnames(coordinates2)<-c("chr","start","end")
    coordinates2[,2]<-as.numeric(coordinates2[,2])
    coordinates2[,3]<-as.numeric(coordinates2[,3])
  } else {
    binaryAccessibility<-ifelse(accessibilityData>0,1,0)
    coordinates2<-read.table(text = pre_scATACseq[,1], sep = " ", colClasses = "character")
    colnames(coordinates2)<-c("chr","start","end")
    coordinates2[,2]<-as.numeric(coordinates2[,2])
    coordinates2[,3]<-as.numeric(coordinates2[,3])
  }
  accessibilityMatrix<-cbind(coordinates2,binaryAccessibility)
  return(accessibilityMatrix)
}