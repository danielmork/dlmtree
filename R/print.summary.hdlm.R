#' print.summary.hdlm
#'
#' @title Prints an overview with summary of model class 'hdlm'
#' @description Method for printing an overview with summary of model class 'hdlm'
#'
#' @param x an object of type 'summary.hdlm', result of call to summary.hdlm()
#' @param digits integer number of digits to round
#' @param cw.only print only results for exposures with critical windows
#' @param ... additional parameters
#'
#' @returns output of hdlm fit in R console
#' @export
#'
print.summary.hdlm <- function(x, digits = 3, cw.only = TRUE, ...)
{
  cat("---\n")
  cat("HDLM summary\n\n")

  # Print model info
  cat("Model run info:\n")
  cat("-", Reduce(paste, deparse1(x$formula)), "\n")
  cat("- sample size:", format(x$n, big.mark = ","), "\n")
  cat("- family:", x$ctr$response, "\n")
  cat("-", x$ctr$n.trees, "trees\n")
  cat("-", x$ctr$n.burn, "burn-in iterations\n")
  cat("-", x$ctr$n.iter, "post-burn iterations\n")
  cat("-", x$ctr$n.thin, "thinning factor\n")
  cat("-", x$modPrior, "modifier sparsity prior\n")
  cat("-", x$conf.level, "confidence level\n")

  # Print fixed effect coefficient results (logistic)
  cat("\nFixed effects:\n")
  if (length(x$droppedCovar) > 0) {
    cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
  }
    
  gamma.out <- data.frame("Mean" = round(x$gamma.mean, digits),
                          "Lower" = round(x$gamma.ci[1,], digits),
                          "Upper" = round(x$gamma.ci[2,], digits))

  row.names(gamma.out) <-
    ifelse(x$gamma.ci[1,] > 0 | x$gamma.ci[2,] < 0,
            paste0("*", names(x$gamma.mean)),
            paste0(" ", names(x$gamma.mean)))
  print(gamma.out)

  cat("---\n")
  cat("* = CI does not contain zero\n")


  cat("\nModifiers:\n")
  print(x$pip)
  cat("---\n")
  cat("PIP = Posterior inclusion probability\n")

  if(x$ctr$response == "gaussian"){
    cat("\nresidual standard errors: ")
    cat(round(x$rse, 3))
  }

  cat("\n---\n")

  cat("To obtain exposure effect estimates, use the 'shiny(fit)' function.\n")

  cat("\n")
}
