#' print.summary.hdlm
#'
#' @param object an object of type 'summary.hdlm', result of call to summary.hdlm()
#' @param digits integer number of digits to round
#' @param cw.only print only results for exposures with critical windows
#'
#' @return output in R console
#' @export print.summary.hdlm
#' @export
#'
print.summary.hdlm <- function(object, digits = 4, cw.only = TRUE)
{
  cat("\nHDLM:\n\n")

  # Print model info
  cat("Model run info:\n")

  # Print ZI and NB part separately for ZINB
  # if (object$family == "zinb") {
  #   cat("- ZI:", Reduce(paste, deparse(object$formula.zi)), "\n")
  #   cat("- NB:", Reduce(paste, deparse(object$formula)), "\n")
  # } else {
  cat("-", Reduce(paste, deparse(object$formula)), "\n")
  # }
  cat("- family:", object$ctr$response, "\n")
  cat("-", object$ctr$n.trees, "trees\n")
  cat("-", object$ctr$n.burn, "burn-in iterations\n")
  cat("-", object$ctr$n.iter, "post-burn iterations\n")
  cat("-", object$ctr$n.thin, "thinning factor\n")
  cat("-", object$modPrior, "modifier sparsity prior\n")
  cat("-", object$conf.level, "confidence level\n")

  # Print fixed effect coefficient results (logistic)
  cat("\nFixed effects:\n")
  if (length(object$droppedCovar) > 0) {
    cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
  }
    
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


  cat("\nModifiers:\n")
  print(object$pip)
  cat("---\n")
  cat("PIP = Posterior inclusion probability\n")

  # # Print fixed effect coefficient results (ZINB - Binary)
  # if (object$family == "zinb") {
  #   cat("\nFixed effects (ZI model):\n")
  #   if (length(object$droppedCovar) > 0) {
  #     cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
  #   }
      
  #   b1.out <- data.frame("Mean" = round(object$b1.mean, digits),
  #                        "Lower" = round(object$b1.ci[1,], digits),
  #                         "Upper" = round(object$b1.ci[2,], digits))

  #   row.names(b1.out) <-
  #     ifelse(object$b1.ci[1,] > 0 | object$b1.ci[2,] < 0,
  #           paste0("*", names(object$b1.mean)),
  #           paste0(" ", names(object$b1.mean)))
  #   print(b1.out)
  #   cat("---\n")

  #   # Print fixed effect coefficient results (ZINB - Count)
  #   cat("\nFixed effects (NB model):\n")
  #   if (length(object$droppedCovar) > 0) {
  #     cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
  #   }
  #   b2.out <- data.frame("Mean" = round(object$b2.mean, digits),
  #                        "Lower" = round(object$b2.ci[1,], digits),
  #                        "Upper" = round(object$b2.ci[2,], digits))
  #   row.names(b2.out) <-
  #     ifelse(object$b2.ci[1,] > 0 | object$b2.ci[2,] < 0,
  #           paste0("*", names(object$b2.mean)),
  #           paste0(" ", names(object$b2.mean)))
  #   print(b2.out)
  #   cat("---\n")
  #   cat("* = CI does not contain zero\n")

  #   # Print dispersion parameter, r
  #   cat("\nFixed effects (Dispersion):\n")
  #   if (length(object$droppedCovar) > 0) {
  #     cat("dropped collinear covariates:", paste(object$droppedCovar, collapse = ", "),"\n")
  #   }
      
  #   r.out <- data.frame("Mean" = round(object$r.mean, digits),
  #                         "Lower" = round(object$r.ci[1], digits),
  #                         "Upper" = round(object$r.ci[2], digits))
  #   row.names(r.out) <- "Dispersion"
  #     #ifelse(object$r.ci[1,] > 0 | object$b2.ci[2,] < 0,
  #     #      paste0("*", names(object$b2.mean)),
  #     #      paste0(" ", names(object$b2.mean)))
  #   print(r.out)
  #   cat("---\n")
  # }


  # Shiny app message here
  # cat("\n--\nExposure effects: critical windows\n")
  # cat("* = Exposure selected by Bayes Factor\n")
  # cat("(x.xx) = Relative effect size\n")
  # # cat("exposure name (signal): critical windows")
  # for (ex.name in names(object$TreeStructs)) {
  #   if (any(object$TreeStructs[[ex.name]]$marg.cw) | !cw.only | object$expSel[ex.name]) {
  #     cat("\n",
  #         paste0(ifelse(object$expSel[ex.name], "*", " "), ex.name,
  #                 " (", round(object$expVar[2, ex.name], 2), "): ",
  #                   ppRange(which(object$TreeStructs[[ex.name]]$marg.cw))))
  #   }
  # }

  cat("\n---\n")

  cat("To obtain exposure effect estimates, use the 'shiny(fit)' function.\n")

  cat("\n")
}
