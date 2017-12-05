#' Distance from Centroid in Principal Coordinate Space using Angular Cosine Distance
#' 
#' When using default parameters, calculates cell-to-cell variation (between columns) in accessibility using angular cosine distance inputted into principal coordinate analysis. The centroid of the cells in principal coordinate space is found, then each cell's distance from the centroid constitutes the variability.
#' @param data An accessibility matrix, with columns as cells and rows as DNA sequences.
#' @export
#' @examples 
#' centroid_distance()

centroid_distance<-function(data, distance = "cosine", method = "centroid") {
  if (distance == "cosine") {
    dist<-cosine_distance(data)} 
  else if (distance == "manhattan") {
    dist<-dist(t(data),method="manhattan")
  }
  dist_centroid<-vegan::betadisper(dist,replicate(ncol(data),1),method)
  return(dist_centroid$distances)
}