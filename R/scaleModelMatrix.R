#' scaleModelMatrix
#'
#' @title Centers and scales a matrix
#' @description Method for centering and scaling a matrix
#'
#' @param M a matrix to center and scale
#'
#' @returns a scaled matrix
#' @export
#'
scaleModelMatrix <- function(M)
{
  if (is.vector(M)) {
    vec <- TRUE
    M   <- as.matrix(M)
  } else {
    vec <- FALSE
  }

  M.center  <- sapply(1:ncol(M), function(j) ifelse(diff(range(M[,j])) > 0 & length(unique(M[,j])) > 2, mean(M[,j]), 0))
  M.scale   <- sapply(1:ncol(M), function(j) sqrt(sum((M[,j] - M.center[j])^2)))
  M         <- scale(M, center = M.center, scale = M.scale)

  if (vec) {
    M <- as.vector(M)
  }

  return(M)
}
