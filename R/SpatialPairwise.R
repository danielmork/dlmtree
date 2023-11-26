#' spatial_pairwise
#' @title Spatial adjacency matrix to two nodes vectors
#' @param ... Adjacency matrix
#'
#' @return Two nodes vectors
#' @export
#'
# function 2: Mirrors the lower triangular to upper triangular aka make is symmetric
spatial_pairwise <- function(AdjMatrix){
  # parameters
  i = 1
  nodes1 = rep(NA, sum(AdjMatrix)/2)
  nodes2 = rep(NA, sum(AdjMatrix)/2)
  
  for(row in 1:nrow(AdjMatrix)){
    edges.index = which(AdjMatrix[row,] == 1)
    if(length(edges.index) != 0){
      for(col in 1:length(edges.index)){
        if(row < edges.index[col]){
          nodes1[i] = row
          nodes2[i] = edges.index[col]
          i = i + 1
        }
      }
    }
  }
  
  return(list("nodes1" = nodes1, 
              "nodes2" = nodes2))
}