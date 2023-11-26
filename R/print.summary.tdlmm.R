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

  # Print ZI and NB part separately for ZINB
  if(object$family == "zinb"){
    cat("- ZI:", Reduce(paste, deparse(object$formula_zi)), "\n")
    cat("- NB:", Reduce(paste, deparse(object$formula)), "\n")
  } else {
    cat("-", Reduce(paste, deparse(object$formula)), "\n")
  }
  cat("- family:", object$family, "\n")
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


  # Print fixed effect coefficient results (logistic)
  if(object$family != "zinb"){
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
  }

  # Print fixed effect coefficient results (ZINB - Binary)
  if(object$family == "zinb"){
    cat("\nFixed effects (Binary):\n")
    if (length(object$droppedCovar) > 0)
      cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
    b1.out <- data.frame("Mean" = round(object$b1.mean, digits),
                         "Lower" = round(object$b1.ci[1,], digits),
                          "Upper" = round(object$b1.ci[2,], digits))
    row.names(b1.out) <-
      ifelse(object$b1.ci[1,] > 0 | object$b1.ci[2,] < 0,
            paste0("*", names(object$b1.mean)),
            paste0(" ", names(object$b1.mean)))
    print(b1.out)
    cat("---\n")

    # Print fixed effect coefficient results (ZINB - Count)
    cat("\nFixed effects (Count):\n")
    if (length(object$droppedCovar) > 0)
      cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
    b2.out <- data.frame("Mean" = round(object$b2.mean, digits),
                          "Lower" = round(object$b2.ci[1,], digits),
                          "Upper" = round(object$b2.ci[2,], digits))
    row.names(b2.out) <-
      ifelse(object$b2.ci[1,] > 0 | object$b2.ci[2,] < 0,
            paste0("*", names(object$b2.mean)),
            paste0(" ", names(object$b2.mean)))
    print(b2.out)
    cat("---\n")
    cat("* = CI does not contain zero\n")

    # Print dispersion parameter, r
    cat("\nFixed effects (Dispersion):\n")
    if (length(object$droppedCovar) > 0)
      cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
    r.out <- data.frame("Mean" = round(object$r.mean, digits),
                          "Lower" = round(object$r.ci[1], digits),
                          "Upper" = round(object$r.ci[2], digits))
    row.names(r.out) <- "Dispersion"
      #ifelse(object$r.ci[1,] > 0 | object$b2.ci[2,] < 0,
      #      paste0("*", names(object$b2.mean)),
      #      paste0(" ", names(object$b2.mean)))
    print(r.out)
    cat("---\n")
  }


  # Print exposure effects
  cat("\n--\nExposure effects: critical windows\n")
  cat("* = Exposure selected by Bayes Factor\n")
  cat("(x.xx) = Relative effect size\n")
  # cat("exposure name (signal): critical windows")
  for (ex.name in names(object$DLM)) {
    if (any(object$DLM[[ex.name]]$marg.cw) | !cw.only | object$expSel[ex.name]) {
      cat("\n",
          paste0(ifelse(object$expSel[ex.name], "*", " "), ex.name,
                 " (", round(object$expVar[2, ex.name], 2), "): ",
                 ppRange(which(object$DLM[[ex.name]]$marg.cw))))
    }
  }


  # Print mixture effects
  if (object$interaction > 0) {
    cat("\n--\nInteraction effects: critical windows\n")
    
    for (mix.name in names(object$MIX)) {
      cw <- rowSums(object$MIX[[mix.name]]$cw)
  
      if (any(cw > 0) | (!cw.only & length(object$MIX) > 1)) {
        if(length(names(object$MIX)) == 1){ # Single interaction when we fit TDLMMns to two components
          cat("\n",
            paste0(ifelse(object$mixSel[mix.name], "*", " "),
                   object$MIX[[mix.name]]$rows, "/", object$MIX[[mix.name]]$cols,
                   " (", round(object$mixVar[1], 2), "):"))
          for (r in which(cw > 0)) {
            s <- which(object$MIX[[mix.name]]$cw[r,])
            cat("\n", paste0(r, "/", ppRange(s)))
          }
        } else { # More than one interaction
          cat("\n",
            paste0(ifelse(object$mixSel[mix.name], "*", " "),
                   object$MIX[[mix.name]]$rows, "/", object$MIX[[mix.name]]$cols,
                   " (", round(object$mixVar[2, mix.name], 2), "):"))
          for (r in which(cw > 0)) {
            s <- which(object$MIX[[mix.name]]$cw[r,])
            cat("\n", paste0(r, "/", ppRange(s)))
          }
        }
        
      } else {
        cat("\n - No critical windows")
      }
    }
  }
  cat("\n---\n")
}
