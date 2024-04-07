#' print.summary.tdlmm
#'
#' @param x an object of type 'summary.tdlmm', result of call to summary.tdlmm()
#' @param digits integer number of digits to round
#' @param cw.only print only results for exposures with critical windows
#' @param ... additional parameters
#'
#' @return output in R console
#' @export print.summary.tdlmm
#' @export
#'
print.summary.tdlmm <- function(x, digits = 3, cw.only = TRUE, ...)
{
  cat("\nTDLMM:\n\n")

  # Print model info
  cat("Model run info:\n")

  # Print ZI and NB part separately for ZINB
  if (x$family == "zinb") {
    cat("- ZI:", Reduce(paste, deparse(x$formula.zi)), "\n")
    cat("- NB:", Reduce(paste, deparse(x$formula)), "\n")
  } else {
    cat("-", Reduce(paste, deparse(x$formula)), "\n")
  }
  cat("- family:", x$family, "\n")
  cat("- ", x$nTrees, " trees (alpha = ", x$treePrior[1], ", beta = ", x$treePrior[2], ")\n", sep = "")
  cat("-", x$nBurn, "burn-in iterations\n")
  cat("-", x$nIter, "post-burn iterations\n")
  cat("-", x$nThin, "thinning factor\n")
  cat("-", x$nExp, "exposures measured at", x$nLags, "time points\n")
  if (x$interaction > 0) {
    cat("-", x$nMix, "two-way interactions")
    if (x$interaction == 1) {
      cat(" (no-self interactions)\n")
    } else {
      cat(" (all interactions)\n")
    }
  }
  cat("-", x$mixPrior, "kappa sparsity prior\n")
  cat("-", x$conf.level, "confidence level\n")


  # Print fixed effect coefficient results (logistic)
  if (x$family != "zinb") {
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
  }

  # Print fixed effect coefficient results (ZINB - Binary)
  if (x$family == "zinb") {
    cat("\nFixed effects (ZI model):\n")
    if (length(x$droppedCovar) > 0) {
      cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
    }
      
    b1.out <- data.frame("Mean" = round(x$b1.mean, digits),
                         "Lower" = round(x$b1.ci[1,], digits),
                          "Upper" = round(x$b1.ci[2,], digits))

    row.names(b1.out) <-
      ifelse(x$b1.ci[1,] > 0 | x$b1.ci[2,] < 0,
            paste0("*", names(x$b1.mean)),
            paste0(" ", names(x$b1.mean)))
    print(b1.out)
    cat("---\n")

    # Print fixed effect coefficient results (ZINB - Count)
    cat("\nFixed effects (NB model):\n")
    if (length(x$droppedCovar) > 0) {
      cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
    }
    b2.out <- data.frame("Mean" = round(x$b2.mean, digits),
                         "Lower" = round(x$b2.ci[1,], digits),
                         "Upper" = round(x$b2.ci[2,], digits))
    row.names(b2.out) <-
      ifelse(x$b2.ci[1,] > 0 | x$b2.ci[2,] < 0,
            paste0("*", names(x$b2.mean)),
            paste0(" ", names(x$b2.mean)))
    print(b2.out)
    cat("---\n")
    cat("* = CI does not contain zero\n")

    # Print dispersion parameter, r
    cat("\nFixed effects (Dispersion):\n")
    if (length(x$droppedCovar) > 0) {
      cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
    }
      
    r.out <- data.frame("Mean" = round(x$r.mean, digits),
                          "Lower" = round(x$r.ci[1], digits),
                          "Upper" = round(x$r.ci[2], digits))
    row.names(r.out) <- "Dispersion"
      #ifelse(x$r.ci[1,] > 0 | x$b2.ci[2,] < 0,
      #      paste0("*", names(x$b2.mean)),
      #      paste0(" ", names(x$b2.mean)))
    print(r.out)
    cat("---\n")
  }


  # Print exposure effects
  cat("\n--\nExposure effects: critical windows\n")
  cat("* = Exposure selected by Bayes Factor\n")
  cat("(x.xx) = Relative effect size\n")
  # cat("exposure name (signal): critical windows")
  for (ex.name in names(x$DLM)) {
    if (any(x$DLM[[ex.name]]$marg.cw) | !cw.only | x$expSel[ex.name]) {
      cat("\n",
          paste0(ifelse(x$expSel[ex.name], "*", " "), ex.name,
                 " (", round(x$expVar[2, ex.name], 2), "): ",
                 ppRange(which(x$DLM[[ex.name]]$marg.cw))))
    }
  }


  # Print mixture effects
  if (x$interaction > 0) {
    cat("\n--\nInteraction effects: critical windows\n")

    cw.any = FALSE # Counting for "no interaction" messages
    for (mix.name in names(x$MIX)) {
      cw <- rowSums(x$MIX[[mix.name]]$cw)
      cw.any <- any(cw > 0)
  
      if (any(cw > 0) | (!cw.only & length(x$MIX) > 1)) {
        if (length(names(x$MIX)) == 1) { # Single interaction when we fit TDLMMns to two components
          cat("\n",
            paste0(ifelse(x$mixSel[mix.name], "*", " "),
                   x$MIX[[mix.name]]$rows, "/", x$MIX[[mix.name]]$cols,
                   " (", round(x$mixVar[1], 2), "):"))
          for (r in which(cw > 0)) {
            s <- which(x$MIX[[mix.name]]$cw[r,])
            cat("\n", paste0(r, "/", ppRange(s)))
          }
        } else { # More than one interaction
          cat("\n",
            paste0(ifelse(x$mixSel[mix.name], "*", " "),
                   x$MIX[[mix.name]]$rows, "/", x$MIX[[mix.name]]$cols,
                   " (", round(x$mixVar[2, mix.name], 2), "):"))
          for (r in which(cw > 0)) {
            s <- which(x$MIX[[mix.name]]$cw[r,])
            cat("\n", paste0(r, "/", ppRange(s)))
          }
        }
        
      } 
    }

    if(!cw.any) {
      cat("\n - No critical windows")
    }
  }

  cat("\n---\n")
}
