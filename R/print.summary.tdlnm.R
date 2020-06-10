#' print.summary.tdlnm
#'
#' @param object an object of type 'summary.tdlnm', result of call to summary.tdlnm()
#' @param digits integer number of digits to round
#'
#' @return
#' @export
#'
#' @examples
print.summary.tdlnm <- function(object, digits = 3)
{
  cat("TDLNM summary\n\n")
  cat("Model run info\n")
  cat("-", object$ctr$n.trees, "trees\n")
  cat("-", object$ctr$n.burn, "burn-in iterations\n")
  cat("-", object$ctr$n.iter, "post-burn iterations\n")
  cat("-", object$ctr$n.thin, "thinning factor\n")
  cat("-", object$conf.level, "confidence level\n")
  
  cat("\nFixed effect coefficients:\n")
  gamma.out <- data.frame("Mean" = round(object$gamma.mean, digits),
                          "Lower Bound" = round(object$gamma.ci[1,], digits),
                          "Upper Bound" = round(object$gamma.ci[2,], digits))
  row.names(gamma.out) <- ifelse(object$gamma.ci[1,] > 0 | object$gamma.ci[2,] < 0,
                                 paste0("*", names(object$gamma.mean)),
                                 names(object$gamma.mean))
  print(gamma.out)
  cat("---\n")
  cat("* = CI does not contain zero\n")
  
  cat("\nDLNM effect:")
  cat("\nrange = ["); cat(round(range(object$matfit), 3), sep = ", "); cat("]")
  cat("\nsignal-to-noise =", round(object$sig.to.noise, digits))
  cat("\ncritical windows: ")
  cw <- which((colSums(object$cilower > 0) + colSums(object$ciupper < 0)) > 0)
  cat(cw, "\n")
}