#' summary.tdlm
#'
#' @param object an object of class 'tdlm', result of a call to tdlnm using
#' zero splits (i.e. a linear effect DLM)
#' @param conf.level confidence level for computation of credible intervals
#'
#' @return
#' @export summary.tdlm
#' @export
#'
summary.tdlm <- function(object,
                         conf.level = 0.95)
{
  Lags <- max(object$DLM$tmax)
  Iter <-  max(object$DLM$Iter)
  ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  dlmest <- dlmEst(as.matrix(object$DLM)[,-c(3:4)], Lags, Iter)

  # DLM Estimates
  matfit <- rowMeans(dlmest)
  cilower <- apply(dlmest, 1, quantile, probs = ci.lims[1])
  ciupper <- apply(dlmest, 1, quantile, probs = ci.lims[2])

  # Cumulative effect estimates
  ce <- colSums(dlmest)
  cumulative.effect <- list("mean" = mean(ce),
                            "ci" = quantile(ce, ci.lims))

  # Fixed effect estimates
  gamma.mean <- colMeans(object$gamma)
  gamma.ci <- apply(object$gamma, 2, quantile, probs = ci.lims)

  # Return
  ret <- list("ctr" = list(dl.function = object$dlFunction,
                           n.trees = object$nTrees,
                           n.iter = object$nIter,
                           n.thin = object$nThin,
                           n.burn = object$nBurn,
                           response = object$family),
              "conf.level" = conf.level,
              "sig.to.noise" = ifelse(is.null(object$sigma2), NA,
                                      var(object$fhat) / mean(object$sigma2)),
              "matfit" = matfit,
              "cilower" = cilower,
              "ciupper" = ciupper,
              "cumulative.effect" = cumulative.effect,
              "gamma.mean" = gamma.mean,
              "gamma.ci" = gamma.ci)
  class(ret) <- "summary.tdlnm"
  return(ret)
}
