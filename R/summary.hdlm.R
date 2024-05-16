#' summary.hdlm
#'
#' @title Creates a summary object of class 'hdlm'
#' @description Method for creating a summary object of class 'hdlm'
#'
#' @param object an object of class 'hdlm'
#' @param conf.level confidence level for computation of credible intervals
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
#' summary(fit)
#'
#' @returns list of type 'summary.hdlm'
#' @export
#'
summary.hdlm <- function(object, conf.level = 0.95, ...)
{
  Lags    <- max(object$TreeStructs$tmax)
  Iter    <- max(object$TreeStructs$Iter)
  ci.lims <-  c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  # Fixed effect estimates
  gamma.mean  <- colMeans(object$gamma)
  gamma.ci    <- apply(object$gamma, 2, quantile, probs = ci.lims)

  # posterior inclusion probability of modifiers
  pip_df <- data.frame("PIP" = pip(object))

  # Return
  ret <- list("ctr" = list(class    = object$class,
                           n.trees  = object$nTrees,
                           n.iter   = object$nIter,
                           n.thin   = object$nThin,
                           n.burn   = object$nBurn,
                           response = object$family),
              "conf.level"   = conf.level,
              "sig.to.noise" = ifelse(is.null(object$sigma2), NA,
                                        var(object$fhat) / mean(object$sigma2)),
              "rse"          = sd(object$sigma2),
              "n"            = nrow(object$data),
              "modPrior"     = object$zeta,
              "pip"          = pip_df,
              "gamma.mean"   = gamma.mean,
              "gamma.ci"     = gamma.ci,
              "formula"      = object$formula)

  class(ret) <- "summary.hdlm"
  
  return(ret)
}
