#' print.summary.tdlnm
#'
#' @param object an object of type 'summary.tdlnm', result of call to summary.tdlnm()
#' @param digits integer number of digits to round
#'
#' @return output in R console
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
  if(object$ctr$response != "zinb"){
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
  }

  # Print fixed effect coefficient results (ZINB - Binary)
  if(object$ctr$response == "zinb"){
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
  cat(cw, "\n")
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
