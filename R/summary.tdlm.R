#' summary.tdlm
#'
#' @title Creates a summary object of class 'tdlm'
#' @description Method for creating a summary object of class 'tdlm'
#'
#' @param object an object of dlm class 'tdlm' (i.e. a linear effect DLM)
#' @param conf.level confidence level for computation of credible intervals
#' @param ... additional parameters
#'
#' @examples
#' D <- sim.tdlmm(sim = "A", mean.p = 0.5, n = 1000)
#' fit <- dlmtree(y ~ .,
#'                data = D$dat,
#'                exposure.data = D$exposures[[1]],
#'                dlm.type = "linear",
#'                family = "logit",
#'                binomial.size = 1)
#' summary(fit)
#'
#' @returns list of type 'summary.tdlm'
#' @export
#'
summary.tdlm <- function(object, conf.level = 0.95, ...){
  Lags    <- max(object$TreeStructs$tmax)
  Iter    <- max(object$TreeStructs$Iter)
  ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  dlmest  <- dlmEst(as.matrix(object$TreeStructs)[,-c(3:4)], Lags, Iter)

  # DLM Estimates
  matfit  <- rowMeans(dlmest)
  cilower <- apply(dlmest, 1, quantile, probs = ci.lims[1])
  ciupper <- apply(dlmest, 1, quantile, probs = ci.lims[2])

  # Cumulative effect estimates
  ce                <- colSums(dlmest)
  cumulative.effect <- c("mean" = mean(ce), quantile(ce, ci.lims))
  xvals             <- seq(object$Xrange[1], object$Xrange[2], length.out = 50)
  cumulative.effect <- data.frame("vals" = xvals,
                                  "mean" = cumulative.effect[1] * xvals,
                                  "lower" = cumulative.effect[2] * xvals,
                                  "upper" = cumulative.effect[3] * xvals)

  # Fixed effect estimates
  gamma.mean  <- colMeans(object$gamma)
  gamma.ci    <- apply(object$gamma, 2, quantile, probs = ci.lims)

  # ZINB
  # binary
  b1.mean     <- colMeans(object$b1)
  b1.ci       <- apply(object$b1, 2, quantile, probs = ci.lims)

  # count
  b2.mean     <- colMeans(object$b2)
  b2.ci       <- apply(object$b2, 2, quantile, probs = ci.lims)

  # Dispersion parameter
  r.mean      <- mean(object$r)
  r.ci        <- quantile(object$r, probs = ci.lims)

  # Return
  ret <- list("ctr" = list(class    = object$class,
                           n.trees  = object$nTrees,
                           n.iter   = object$nIter,
                           n.thin   = object$nThin,
                           n.burn   = object$nBurn,
                           response = object$family),
              "conf.level"        = conf.level,
              "sig.to.noise"      = ifelse(is.null(object$sigma2), NA,
                                        var(object$fhat) / mean(object$sigma2)),
              "matfit"            = matfit,
              "cilower"           = cilower,
              "ciupper"           = ciupper,
              "cumulative.effect" = cumulative.effect,
              "gamma.mean"        = gamma.mean,
              "gamma.ci"          = gamma.ci,
              "b1.mean"           = b1.mean,
              "b1.ci"             = b1.ci,
              "b2.mean"           = b2.mean,
              "b2.ci"             = b2.ci,
              "r.mean"            = r.mean,
              "r.ci"              = r.ci,
              "formula"           = object$formula,
              "formula.zi"        = object$formula.zi)

  class(ret) <- "summary.tdlm"
  
  return(ret)
}
