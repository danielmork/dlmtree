#' @method summary tdlm
#' @rdname summary
#'
#' @export
summary.tdlm <- function(x, conf.level = 0.95, mcmc = FALSE, ...){
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
              "n.lag"             = Lags,
              "sig.to.noise"      = ifelse(is.null(x$sigma2), NA,
                                        var(x$fhat) / mean(x$sigma2)),
              "rse"               = sd(x$sigma2),
              "n"                 = nrow(x$data),
              "matfit"            = matfit,
              "cilower"           = cilower,
              "ciupper"           = ciupper,
              "cumulative.effect" = cumulative.effect,
              "formula"           = x$formula)
  

  if (!x$zinb) {
    # Gaussian / Logistic
    ret$gamma.mean  <- gamma.mean
    ret$gamma.ci    <- gamma.ci

  } else {
    # ZINB
    ret$formula.zi  <- x$formula.zi
    
    # binary
    ret$b1.mean <- b1.mean
    ret$b1.ci   <- b1.ci
    
    # count
    ret$b2.mean <- b2.mean
    ret$b2.ci   <- b2.ci
    
    # Dispersion parameter
    ret$r.mean  <- r.mean
    ret$r.ci    <- r.ci
  }

  
  mcmc.samples <- list()    
  if (mcmc) {
    # dlm
    mcmc.samples$dlm.mcmc           <- t(dlmest)
    colnames(mcmc.samples$dlm.mcmc) <- 1:Lags
    mcmc.samples$cumulative.mcmc    <- data.frame("cumulative effects" = ce)
    
    # tree
    mcmc.samples$tree.size <- x$termNodes
    
    # mhr parameters
    mcmc.samples$accept <- x$treeAccept[, 1:2]
    colnames(mcmc.samples$accept) <- c("step", "success")
    
    
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

  class(ret) <- "summary.tdlm"
  
  return(ret)
}


