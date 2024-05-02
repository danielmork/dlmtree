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
#' @examples
#' D <- sim.hdlmm(sim = "B", n = 1000)
#' fit <- dlmtree(y ~ ., 
#'                data = D$dat,
#'                exposure.data = D$exposures,
#'                dlm.type = "linear",
#'                family = "gaussian",
#'                het = TRUE)
#' fit_sum <- summary(fit)
#' print(fit_sum)
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

  # Print ZI and NB part separately for ZINB
  # if (x$family == "zinb") {
  #   cat("- ZI:", Reduce(paste, deparse(x$formula.zi)), "\n")
  #   cat("- NB:", Reduce(paste, deparse(x$formula)), "\n")
  # } else {
  cat("-", Reduce(paste, deparse(x$formula)), "\n")
  # }
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

  # # Print fixed effect coefficient results (ZINB - Binary)
  # if (x$family == "zinb") {
  #   cat("\nFixed effects (ZI model):\n")
  #   if (length(x$droppedCovar) > 0) {
  #     cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
  #   }
      
  #   b1.out <- data.frame("Mean" = round(x$b1.mean, digits),
  #                        "Lower" = round(x$b1.ci[1,], digits),
  #                         "Upper" = round(x$b1.ci[2,], digits))

  #   row.names(b1.out) <-
  #     ifelse(x$b1.ci[1,] > 0 | x$b1.ci[2,] < 0,
  #           paste0("*", names(x$b1.mean)),
  #           paste0(" ", names(x$b1.mean)))
  #   print(b1.out)
  #   cat("---\n")

  #   # Print fixed effect coefficient results (ZINB - Count)
  #   cat("\nFixed effects (NB model):\n")
  #   if (length(x$droppedCovar) > 0) {
  #     cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
  #   }
  #   b2.out <- data.frame("Mean" = round(x$b2.mean, digits),
  #                        "Lower" = round(x$b2.ci[1,], digits),
  #                        "Upper" = round(x$b2.ci[2,], digits))
  #   row.names(b2.out) <-
  #     ifelse(x$b2.ci[1,] > 0 | x$b2.ci[2,] < 0,
  #           paste0("*", names(x$b2.mean)),
  #           paste0(" ", names(x$b2.mean)))
  #   print(b2.out)
  #   cat("---\n")
  #   cat("* = CI does not contain zero\n")

  #   # Print dispersion parameter, r
  #   cat("\nFixed effects (Dispersion):\n")
  #   if (length(x$droppedCovar) > 0) {
  #     cat("dropped collinear covariates:", paste(x$droppedCovar, collapse = ", "),"\n")
  #   }
      
  #   r.out <- data.frame("Mean" = round(x$r.mean, digits),
  #                         "Lower" = round(x$r.ci[1], digits),
  #                         "Upper" = round(x$r.ci[2], digits))
  #   row.names(r.out) <- "Dispersion"
  #     #ifelse(x$r.ci[1,] > 0 | x$b2.ci[2,] < 0,
  #     #      paste0("*", names(x$b2.mean)),
  #     #      paste0(" ", names(x$b2.mean)))
  #   print(r.out)
  #   cat("---\n")
  # }


  # Shiny app message here
  # cat("\n--\nExposure effects: critical windows\n")
  # cat("* = Exposure selected by Bayes Factor\n")
  # cat("(x.xx) = Relative effect size\n")
  # # cat("exposure name (signal): critical windows")
  # for (ex.name in names(x$TreeStructs)) {
  #   if (any(x$TreeStructs[[ex.name]]$marg.cw) | !cw.only | x$expSel[ex.name]) {
  #     cat("\n",
  #         paste0(ifelse(x$expSel[ex.name], "*", " "), ex.name,
  #                 " (", round(x$expVar[2, ex.name], 2), "): ",
  #                   ppRange(which(x$TreeStructs[[ex.name]]$marg.cw))))
  #   }
  # }

  cat("\n---\n")

  cat("To obtain exposure effect estimates, use the 'shiny(fit)' function.\n")

  cat("\n")
}
