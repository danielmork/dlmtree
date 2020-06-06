#' Scale model matrix
#'
#' @param M Matrix of predictors
#'
#' @return Matrix of scaled predictors
#'
#' @examples
scaleModelMatrix <- function(M)
{
  if (is.vector(M)) {
    vec <- T
    M <- as.matrix(M)
  } else {
    vec <- F
  }

  # M.scale <- sapply(1:ncol(M), function(j) ifelse(diff(range(M[,j])) > 0, diff(range(M[,j])), 1))
  M.scale <- sapply(1:ncol(M), function(j) sqrt(crossprod(M[,j])))
  M.center <- sapply(1:ncol(M), function(j) ifelse(diff(range(M[,j])) > 0, mean(M[,j]), 0))#ifelse(length(unique(M[,j])) > 2, mean(M[,j]), 0))
  M <- scale(M, M.center, M.scale)

  if (vec) {
    M.scale <- as.numeric(M.scale)
    M.center <- as.numeric(M.center)
    M <- as.vector(M)
  }

  return(M)
}
