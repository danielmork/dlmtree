#' print.summary.tdlnm
#'
#' @param x an object of type 'summary.tdlnm', result of call to summary.tdlnm()
#' @param digits integer number of digits to round
#' @param ... additional parameters
#'
#' @examples
#' D <- sim.tdlnm(sim = "A", error.to.signal = 1)
#' fit <- dlmtree(formula = y ~ .,
#'                data = D$dat,
#'                exposure.data = as.matrix(D$exposures),
#'                dlm.type = "nonlinear",
#'                family = "gaussian")
#' fit_sum <- summary(fit)
#' print(fit_sum)
#' 
#' @returns output of tdlnm fit in R console
#' @export
#'
print.summary.tdlnm <- function(x, digits = 3, ...)
{
  cat("---\n")
  cat("TDLNM summary\n\n")

  cat("Model run info:\n")
  # Print ZI and NB part separately for ZINB
  if (x$ctr$response == "zinb") {
    cat("- ZI:", Reduce(paste, deparse(x$formula.zi)), "\n")
    cat("- NB:", Reduce(paste, deparse(x$formula)), "\n")
  } else {
    cat("-", Reduce(paste, deparse(x$formula)), "\n")
  }

  cat("- family:", x$ctr$response, "\n")
  cat("-", x$ctr$n.trees, "trees\n")
  cat("-", x$ctr$n.burn, "burn-in iterations\n")
  cat("-", x$ctr$n.iter, "post-burn iterations\n")
  cat("-", x$ctr$n.thin, "thinning factor\n")
  cat("-", x$conf.level, "confidence level\n")

  if (x$ctr$response != "zinb") {
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
  } else {  # zinb
    # Print fixed effect coefficient results (ZINB - Binary)
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
  

  cat("\nDLNM effect:")
  cat("\nrange = ["); cat(round(range(x$matfit), 3), sep = ", "); cat("]")

  if (!is.na(x$sig.to.noise)) {
    cat("\nsignal-to-noise =", round(x$sig.to.noise, digits))
  }

  cat("\ncritical windows: ")
  cw <- ppRange(which((colSums(x$cilower > 0) + colSums(x$ciupper < 0)) > 0))

  cat(cw, "\n")
}
