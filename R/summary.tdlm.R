#' @method summary tdlm
#' @rdname summary
#'
#' @export
summary.tdlm <- function(x, conf.level = 0.95, ...){
  Lags    <- max(x$TreeStructs$tmax)
  Iter    <- max(x$TreeStructs$Iter)
  ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  dlmest  <- dlmEst(as.matrix(x$TreeStructs)[,-c(3:4)], Lags, Iter)

  # DLM Estimates
  matfit  <- rowMeans(dlmest)
  cilower <- apply(dlmest, 1, quantile, probs = ci.lims[1])
  ciupper <- apply(dlmest, 1, quantile, probs = ci.lims[2])

  # Cumulative effect estimates
  ce                <- colSums(dlmest)
  cumulative.effect <- c("mean" = mean(ce), quantile(ce, ci.lims))
  xvals             <- seq(x$Xrange[1], x$Xrange[2], length.out = 50)
  

  # Fixed effect estimates
  gamma.mean  <- colMeans(x$gamma)
  gamma.ci    <- apply(x$gamma, 2, quantile, probs = ci.lims)

  # ZINB
  # binary
  b1.mean     <- colMeans(x$b1)
  b1.ci       <- apply(x$b1, 2, quantile, probs = ci.lims)

  # count
  b2.mean     <- colMeans(x$b2)
  b2.ci       <- apply(x$b2, 2, quantile, probs = ci.lims)

  # Dispersion parameter
  r.mean      <- mean(x$r)
  r.ci        <- quantile(x$r, probs = ci.lims)

  # Return
  ret <- list("ctr" = list(class    = x$class,
                           n.trees  = x$nTrees,
                           n.iter   = x$nIter,
                           n.thin   = x$nThin,
                           n.burn   = x$nBurn,
                           response = x$family),
              "conf.level"        = conf.level,
              "sig.to.noise"      = ifelse(is.null(x$sigma2), NA,
                                        var(x$fhat) / mean(x$sigma2)),
              "rse"               = sd(x$sigma2),
              "n"                 = nrow(x$data),
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
              "formula"           = x$formula,
              "formula.zi"        = x$formula.zi)

  class(ret) <- "summary.tdlm"
  
  return(ret)
}
