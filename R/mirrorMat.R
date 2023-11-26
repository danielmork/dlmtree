#' mirrorMat
#' @title mirroring the matrix
#' @param ... NA
#'
#' @return A symmetric matrix
#' @export
#'
# function 2: Mirrors the lower triangular to upper triangular aka make is symmetric
mirrorMat <- function(matrix) {
    matrix[upper.tri(matrix)] <- t(matrix)[upper.tri(matrix)]
    return(matrix)
}
