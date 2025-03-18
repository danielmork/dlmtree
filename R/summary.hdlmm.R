#' @method summary hdlmm
#' @rdname summary
#'
#' @export
summary.hdlmm <- function(x, conf.level = 0.95, ...)
{
  nLags   <- max(x$TreeStructs$tmax)
  Iter    <- max(x$TreeStructs$Iter)
  ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  # Fixed effect estimates
  gamma.mean  <- colMeans(x$gamma)
  gamma.ci    <- apply(x$gamma, 2, quantile, probs = ci.lims)

  # posterior inclusion probability of modifiers
  pip_df <- data.frame("PIP" = pip(x))

  # Return
  ret <- list("ctr" = list(class    = x$class,
                           n.trees  = x$nTrees,
                           n.iter   = x$nIter,
                           n.thin   = x$nThin,
                           n.burn   = x$nBurn,
                           response = x$family),
              "nLags"         = nLags,
              "nExp"          = x$nExp,
              "nMix"          = x$nMix,
              "interaction"   = x$interaction,
              "mixPrior"      = x$mixPrior,
              "modPrior"      = x$zeta,
              "conf.level"    = conf.level,
              "sig.to.noise"  = ifelse(is.null(x$sigma2), NA,
                                        var(x$fhat) / mean(x$sigma2)),
              "rse"          = sd(x$sigma2),
              "n"            = nrow(x$data),
              "pip"           = pip_df,
              "gamma.mean"    = gamma.mean,
              "gamma.ci"      = gamma.ci,
              "formula"       = x$formula)

  class(ret) <- "summary.hdlmm"
  
  return(ret)
}
