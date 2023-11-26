#' lowerAdj
#' @title A lower triangle adjacency matrix
#' @param ... NA
#'
#' @return a lower adjacency matrix
#' @export
#'
lowerAdj <- function(v, nEdges=1) {
    edges.max <- v*(v-1)/2 # v choose 2
    # Assert length(v)==1 && 1 <= v
    # Assert 0 <= nEdges <= edges.max
    index.edges <- lapply(list(1:(v-1)), function(k) rep(k*(k+1)/2, v-k)) 
    index.edges <- index.edges[[1]] + 1:edges.max
    graph.adjacency <- matrix(0, ncol=v, nrow=v)
    graph.adjacency[sample(index.edges, nEdges)] <- 1
    return(graph.adjacency)
}