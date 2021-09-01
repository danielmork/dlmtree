#' print.summary.tdlnm
#'
#' @param object an object of type 'summary.tdlnm', result of call to summary.tdlnm()
#' @param digits integer number of digits to round
#'
#' @return
#' @export print.summary.tdlnm
#' @export
#'
print.summary.tdlnm <- function(object, digits = 3)
{
  cat("---\n")
  if (object$ctr$dl.function == "tdlnm")
    cat("TDLNM summary\n\n")
  else
    cat("TDLM summary\n\n")
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

  if (object$ctr$dl.function == "tdlnm")
    cat("\nDLNM effect:")
  else
    cat("\nDLM effect:")
  cat("\nrange = ["); cat(round(range(object$matfit), 3), sep = ", "); cat("]")
  if (!is.na(object$sig.to.noise))
    cat("\nsignal-to-noise =", round(object$sig.to.noise, digits))
  cat("\ncritical windows: ")
  if (object$ctr$dl.function == "tdlnm")
    cw <- which((colSums(object$cilower > 0) + colSums(object$ciupper < 0)) > 0)
  else
    cw <- which(object$cilower > 0 | object$ciupper < 0)
  cat(ppRange(cw), "\n")
  if (object$ctr$dl.function == "tdlm") {
    dlm.out <- data.frame("Mean" = round(object$matfit, digits),
                          "Lower Bound" = round(object$cilower, digits),
                          "Upper Bound" = round(object$ciupper, digits))
    rownames(dlm.out) <- ifelse(object$cilower > 0 | object$ciupper < 0,
                                paste0("*Period ", 1:nrow(dlm.out)),
                                paste0("Period ", 1:nrow(dlm.out)))
    print(dlm.out)
    cat("---\n")
    cat("* = CI does not contain zero\n")
  }
}
