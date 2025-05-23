#' @method summary tdlnm
#' @rdname summary
#'
#' @export
summary.tdlnm <- function(x, conf.level = 0.95, pred.at = NULL, cenval = 0, exposure.se = NULL, mcmc = FALSE, verbose = TRUE, ...)
{
  Iter    <- x$mcmcIter
  Lags    <- x$pExp
  ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  if (is.null(exposure.se) && !is.na(x$SE[1])) {
    exposure.se <- mean(as.matrix(x$SE))
  } else if (is.null(exposure.se)) {
    exposure.se <- 0
  }

  # Determine how to break up data for summary
  if (is.null(pred.at)) {
    pred.at <- x$Xsplits
  } else {
    pred.at <- sort(unique(pred.at[which(pred.at >= x$Xrange[1] & pred.at <= x$Xrange[2])]))
  }
  edge.vals <- sort(c(x$Xrange, rowMeans(cbind(pred.at[-1], pred.at[-length(pred.at)]))))


  # Calculate DLNM estimates for gridded values 'pred.at'
  if (verbose) {
    cat("Centered DLNM at exposure value", cenval, "\n")
  }
  cen.quant <- which.min(abs(pred.at - cenval))

  if (exposure.se == 0) {
    dlmest <- dlnmEst(as.matrix(x$TreeStructs), pred.at, Lags, Iter, cen.quant, exposure.se)
  } else {
    dlmest <- dlnmEst(as.matrix(x$TreeStructs), pred.at, Lags, Iter, cenval, exposure.se)
  }

  # Bayes factor
  splitIter <- matrix(1, Lags, Iter)
  splitProb <- rep(0, Lags)

  # Generate cumulative estimtes
  cumexp <- as.data.frame(t(sapply(1:length(pred.at), function(i) {
    cs <- colSums(dlmest[,i,])
    c(pred.at[i], "mean" = mean(cs), quantile(cs, ci.lims))
  })))
  colnames(cumexp) <- c("vals", "mean", "lower", "upper")

  # Matrix of DLNM surface means and CIs, and plot data
  plot.dat  <- as.data.frame(matrix(0, (Lags * length(pred.at)), 10))
  matfit    <- matrix(0, length(pred.at), Lags)
  cilower   <- matrix(0, length(pred.at), Lags)
  ciupper   <- matrix(0, length(pred.at), Lags)

  rownames(matfit) <- rownames(cilower) <- rownames(ciupper) <- pred.at
  colnames(plot.dat) <-  c("Tmin", "Tmax", "Xmin", "Xmax", "PredVal", "Est", "SD", "CIMin", "CIMax", "Effect")
  for (i in 1:Lags) {
    for (j in 1:length(pred.at)) {
      coordest  <- dlmest[i,j,]  
      me        <- mean(coordest)
      s         <- sd(coordest)
      ci        <- quantile(coordest, ci.lims)
      effect    <- ifelse(min(ci) > 0, 1, ifelse(max(ci) < 0, -1, 0))
      plot.dat[(i - 1) * length(pred.at) + j, ] <-
        c(i - 1, i, edge.vals[j], edge.vals[j + 1], pred.at[j],
          me, s, ci, effect)

      matfit[j, i]  <- mean(coordest)
      cilower[j, i] <- ci[1]
      ciupper[j, i] <- ci[2]
    }
  }
  plot.dat$Effect <- factor(plot.dat$Effect, levels = c(-1, 0, 1), labels = c("-", " ", "+"))

  # Fixed effect estimates
  gamma.mean  <- colMeans(x$gamma)
  gamma.ci    <- apply(x$gamma, 2, quantile, probs = ci.lims)

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
              "plot.dat"          = plot.dat,
              "matfit"            = matfit,
              "cilower"           = cilower,
              "ciupper"           = ciupper,
              "cenval"            = cenval,
              "cumulative.effect" = cumexp,
              "pred.at"           = pred.at,
              "gamma.mean"        = gamma.mean,
              "gamma.ci"          = gamma.ci,
              "splitProb"         = splitProb,
              "splitIter"         = splitIter,
              "formula"           = x$formula)
  
  mcmc.samples <- list()    
  if (mcmc) {
    # dlm
    mcmc.samples$dlm.mcmc                <- dlmest
    dimnames(mcmc.samples$dlm.mcmc)[[1]] <- 1:Lags
    dimnames(mcmc.samples$dlm.mcmc)[[2]] <- pred.at
    
    mcmc.samples$cumulative.mcmc           <- sapply(1:length(pred.at), function(i) { colSums(dlmest[,i,]) })
    colnames(mcmc.samples$cumulative.mcmc) <- as.character(pred.at)
    
    # tree
    mcmc.samples$tree.size <- x$termNodes
    
    # mhr parameters
    mcmc.samples$accept <- x$treeAccept[, 1:2]
    names(mcmc.samples$accept) <- c("step", "success")
    
    
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
  
  class(ret) <- "summary.tdlnm"
  
  return(ret)
}
