#' @method summary hdlm
#' @rdname summary
#'
#' @export
summary.hdlm <- function(x, conf.level = 0.95, mcmc = FALSE, ...)
{
  Lags    <- max(x$TreeStructs$tmax)
  Iter    <- max(x$TreeStructs$Iter)
  ci.lims <-  c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

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
              "conf.level"   = conf.level,
              "n.lag"        = Lags,
              "sig.to.noise" = ifelse(is.null(x$sigma2), NA,
                                        var(x$fhat) / mean(x$sigma2)),
              "rse"          = sd(x$sigma2),
              "n"            = nrow(x$data),
              "modPrior"     = x$zeta,
              "mod.names"    = x$modNames,
              "pip"          = pip_df,
              "gamma.mean"   = gamma.mean,
              "gamma.ci"     = gamma.ci,
              "formula"      = x$formula)

  class(ret) <- "summary.hdlm"
  
  mcmc.samples <- list()    
  if (mcmc) {
    # Used for diagnose function
    ret$modIsNum <- x$modIsNum
    ret$data <- x$data
    
    # dlm
    mcmc.samples$dlm.mcmc     <- x$TreeStructs
    
    # tree
    mcmc.samples$modtree.size <- x$termNodesDLM
    mcmc.samples$dlmtree.size <- x$termNodesMod
    
    # mhr parameters
    mcmc.samples$mod.accept <- x$treeModAccept[, 1:2]
    mcmc.samples$dlm.accept <- x$treeDLMAccept[, 1:2]
    
    colnames(mcmc.samples$mod.accept) <- c("step", "success")
    colnames(mcmc.samples$dlm.accept) <- c("step", "success")
    
    
    # other parameters
    if (x$zinb) {
      param.vec <- c("tau", "nu", "sigma2", "b1", "b2")
    } else {
      param.vec <- c("tau", "nu", "sigma2", "gamma")
    }
    
    hyper <- list()
    for (param in param.vec) {
      if (param %in% names(x)) {
        hyper[[param]] <- x[[param]]
      }
    }
    mcmc.samples$hyper <- do.call(cbind, hyper)
    
    ret$mcmc.samples <- mcmc.samples
  }
  
  return(ret)
}
