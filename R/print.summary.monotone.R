#' print.summary.monotone
#'
#' @title Prints an overview with summary of model class 'monotone'
#' @description Method for printing an overview with summary of model class 'monotone'
#' 
#' @param x an object of type 'summary.monotone', result of call to summary.monotone()
#' @param digits integer number of digits to round
#' @param ... additional parameters
#'
#' @returns output in R console
#' @export
#'
print.summary.monotone <- function(x, digits = 3, ...)
{
  cat("---\n")
  cat("monotone-TDLNM summary\n\n")

  cat("Model run info:\n")
  cat("-", Reduce(paste, deparse1(x$formula)), "\n")

  cat("- sample size:", format(x$n, big.mark = ","), "\n")
  cat("- family:", x$ctr$response, "\n")
  cat("-", x$ctr$n.trees, "trees\n")
  cat("-", x$ctr$n.burn, "burn-in iterations\n")
  cat("-", x$ctr$n.iter, "post-burn iterations\n")
  cat("-", x$ctr$n.thin, "thinning factor\n")
  cat("-", x$conf.level, "confidence level\n")

  cat("\nFixed effect coefficients:\n")
  gamma.out <- data.frame("Mean" = round(x$gamma.mean, digits),
                          "Lower" = round(x$gamma.ci[1,], digits),
                          "Upper" = round(x$gamma.ci[2,], digits))
  row.names(gamma.out) <- ifelse(x$gamma.ci[1,] > 0 | x$gamma.ci[2,] < 0,
                                paste0("*", names(x$gamma.mean)),
                                names(x$gamma.mean))
  print(gamma.out)
  cat("---\n")
  cat("* = CI does not contain zero\n")
  
  cat("\nDLNM effect:")
  cat("\nrange = ["); cat(round(range(x$matfit), 3), sep = ", "); cat("]")

  if (!is.na(x$sig.to.noise)) {
    cat("\nsignal-to-noise =", round(x$sig.to.noise, digits))
  }

  cat("\ncritical windows: ")
  cw <- ppRange(which(x$splitProb >= x$conf.level))
  
  if(length(cw) == 0){
    cat("No critical windows \n")
  } else {
    cat(cw, "\n")
  }

  if(x$ctr$response == "gaussian"){
    cat("\nresidual standard errors: ")
    cat(round(x$rse, 3), "\n")
  }

}