#' logPSplit
#'
#' @param alpha Parameter for splitting probability
#' @param beta Parameter for splitting probability
#' @param node.depth Parameter for splitting probability
#' @param terminal Is node terminal
#'
#' @return log probability of split, or log probability of terminal
#'
#' @examples
logPSplit <- function(alpha, beta, node.depth, terminal = F)#, eq = "geom"
{
  # if (eq == "geom") {
    p <- alpha / ((1 + node.depth)^beta)
    if (terminal) {
      return(log1p(-p))
    } else {
      return(log(p))
    }
  # } else if (eq == "logit") {
  #   e <- 1 + (4 ^ (node.depth - alpha))
  #   if (terminal) {
  #     return(log1p(-(1 / e)))
  #   } else {
  #     return(-log(e))
  #   }
  # }
}
