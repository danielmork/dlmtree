#' @method print summary.hdlmm
#' @rdname print.summary
print.summary.hdlmm <- function(x, digits = 3, ...)
{
  cat("---\n")
  cat("HDLMM summary\n\n")

  # Print model info
  cat("Model run info:\n")
  cat("-", Reduce(paste, deparse1(x$formula)), "\n")
  cat("- family:", x$ctr$response, "\n")
  cat("-", x$ctr$n.trees, "trees\n")
  cat("-", x$ctr$n.burn, "burn-in iterations\n")
  cat("-", x$ctr$n.iter, "post-burn iterations\n")
  cat("-", x$ctr$n.thin, "thinning factor\n")
  cat("-", x$n.exp, "exposures measured at", x$n.lag, "time points\n")
  if (x$interaction > 0) {
    cat("-", x$n.mix, "two-way interactions")
    if (x$interaction == 1) {
      cat(" (no-self interactions)\n")
    } else {
      cat(" (all interactions)\n")
    }
  }

  cat("-", x$mod.prior, "modifier sparsity prior\n")
  cat("-", x$mix.prior, "exposure sparsity prior\n")
  cat("-", x$conf.level, "confidence level\n")

  # Print fixed effect coefficient results (logistic)
  cat("\nFixed effects:\n")
  if (length(x$dropped.covar) > 0) {
    cat("dropped collinear covariates:", paste(x$dropped.covar, collapse = ", "),"\n")
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
