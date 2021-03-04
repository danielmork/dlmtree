#' print.summary.tdlmm
#'
#' @param object an object of type 'summary.tdlmm', result of call to summary.tdlmm()
#' @param digits integer number of digits to round
#' @param cw.only print only results for exposures with critical windows
#'
#' @return output in R console
#' @export print.summary.tdlmm
#' @export
#'
print.summary.tdlmm <- function(object, digits = 4, cw.only = TRUE)
{
  cat("\nTDLMM:\n\n")

  # Print model info
  cat("Model run info\n")
  cat("-", Reduce(paste, deparse(object$formula)), "\n")
  cat("-", object$nTrees, "trees (alpha =", object$treePrior[1], ", beta =", object$treePrior[2], ")\n")
  cat("-", object$nBurn, "burn-in iterations\n")
  cat("-", object$nIter, "post-burn iterations\n")
  cat("-", object$nThin, "thinning factor\n")
  cat("-", object$nExp, "exposures measured at", object$nLags, "time points\n")
  if (object$interaction > 0) {
    cat("-", object$nMix, "two-way interactions")
    if (object$interaction == 1) {
      cat(" (no-self interactions)\n")
    } else {
      cat(" (all interactions)\n")
    }
  }
  cat("-", object$mixPrior, "kappa sparsity prior\n")
  cat("-", object$conf.level, "confidence level\n")


  # Print fixed effect coefficient results
  cat("\nFixed effects:\n")
  if (length(object$droppedCovar) > 0)
    cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
  gamma.out <- data.frame("Mean" = round(object$gamma.mean, digits),
                          "Lower" = round(object$gamma.ci[1,], digits),
                          "Upper" = round(object$gamma.ci[2,], digits))
  row.names(gamma.out) <-
    ifelse(object$gamma.ci[1,] > 0 | object$gamma.ci[2,] < 0,
           paste0("*", names(object$gamma.mean)),
           paste0(" ", names(object$gamma.mean)))
  print(gamma.out)
  cat("---\n")
  cat("* = CI does not contain zero\n")


  # Print exposure effects
  cat("\n--\nExposure effects: critical windows\n")
  cat("* = Exposure selected by Bayes Factor\n")
  cat("(x.xx) = Relative effect size\n")
  # cat("exposure name (signal): critical windows")
  for (ex.name in object$expNames) {
    if (length(which(object$DLM[[ex.name]]$marg.cw)) > 0 | !cw.only |
        object$expSel[ex.name]) {
      cat("\n",
          paste0(ifelse(object$expSel[ex.name], "*", " "), ex.name,
                 " (", round(object$expVar[2, ex.name], 2), "): ",
                 ppRange(which(object$DLM[[ex.name]]$marg.cw))))
    }
  }


  # Print mixture effects
  if (object$interaction > 0) {
    cat("\n--\nInteraction effects: critical windows\n")
    # cat(" exposure1 / exposure2 (signal): critical windows")
    for (mix.name in names(object$MIX)) {
      cw <- rowSums(object$MIX[[mix.name]]$cw)
      if (any(cw > 0) | !cw.only) {
        cat("\n",
            paste0(ifelse(object$mixSel[mix.name], "*", " "),
                   object$MIX[[mix.name]]$rows, "/", object$MIX[[mix.name]]$cols,
                   " (", round(object$mixVar[2, mix.name], 2), "):"))
        for (r in which(cw > 0)) {
          s <- which(object$MIX[[mix.name]]$cw[r,])
          cat("\n", paste0(r, "/", ppRange(s)))
        }
      }
    }
  }
  cat("\n---\n")
}
