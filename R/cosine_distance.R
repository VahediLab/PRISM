#' Angular cosine distance
#' 
#' Calculates angular cosine distance between columns (cells) of an accessibility matrix.
#' @param data An accessibility matrix, with columns as cells and rows as DNA sequences.
#' @examples 
#' cosine_distance()

cosine_distance<-function(data){
  v2<-sweep(data,2,apply(data,2,function(a){sqrt(sum(a^2))}),FUN="/")
  v2<-as.matrix(v2)
  cosalpha<-t(v2)%*%v2
  cosalpha[which(cosalpha>1)]=1
  alpha2<-acos(cosalpha)/(pi/2)
  fin <- diag(ncol(data))
  fin[lower.tri(fin, diag=FALSE)]<-alpha2[upper.tri(alpha2)]
  out<-as.dist(fin)
  out[is.nan(out)] <- 1
  return(out)
}