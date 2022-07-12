#' summary.tdlnm
#'
#' @param object an object of class 'tdlnm', result of a call to tdlnm
#' @param pred.at numerical vector of exposure values to make predictions for
#' at each time period
#' @param cenval scalar exposure value that acts as a reference point for
#' predictions at all other exposure values
#' @param conf.level confidence level for computation of credible intervals
#' @param exposure.se scalar smoothing factor, if different from model
#'
#' @return
#' @export summary.tdlnm
#' @export
#'
summary.tdlnm <- function(object,
                          pred.at = NULL,
                          cenval = 0,
                          conf.level = 0.95,
                          exposure.se = NULL,
                          zirt.cond = FALSE,
                          verbose = TRUE)
{
  Iter <- object$mcmcIter
  Lags <- object$pExp
  if (object$monotone) {
    ci.lims <- c((1 - conf.level), 1)
  } else {
    ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  }
  if (is.null(exposure.se) && !is.na(object$SE[1]))
    exposure.se <- mean(as.matrix(object$SE))
  else if (is.null(exposure.se))
    exposure.se = 0

  # Determine how to break up data for summary
  if (is.null(pred.at)) {
    pred.at <- object$Xsplits
  } else {
    pred.at <- sort(unique(pred.at[which(pred.at >= object$Xrange[1] & pred.at <= object$Xrange[2])]))
  }
  edge.vals <- sort(c(object$Xrange, rowMeans(cbind(pred.at[-1], pred.at[-length(pred.at)]))))



  # Calculate DLNM estimates for gridded values 'pred.at'
  if (verbose)
    cat("Centered DLNM at exposure value", cenval, "\n")
  cen.quant <- which.min(abs(pred.at - cenval))
  if (object$shape == "Piecewise Linear") {
    dlmest <- dlnmPLEst(as.matrix(object$DLM), pred.at, Lags, Iter, cen.quant)
  } else if (exposure.se == 0) {
    dlmest <- dlnmEst(as.matrix(object$DLM), pred.at, Lags, Iter,
                      cen.quant, exposure.se)
  } else {
    dlmest <- dlnmEst(as.matrix(object$DLM), pred.at, Lags, Iter,
                      cenval, exposure.se)
  }

  # Bayes factor
  splitIter <- matrix(1, Lags, Iter)
  splitProb <- logBF <- rep(0, Lags)
  if (object$monotone) {
    # splitIter <- object$timeCounts#
    splitIter <- splitPIP(as.matrix(object$DLM[object$DLM$est > 0,]), Lags, Iter)
    splitProb <- apply(splitIter, 1, mean)
    logBF <- log10(splitProb) - log10(1 - splitProb) -
      (log10(1 - (1 - object$p_t)^object$nTrees) - log10((1 - object$p_t)^object$nTrees))
  }

  # Generate cumulative estimtes
  cumexp <- as.data.frame(t(sapply(1:length(pred.at), function(i) {
    cs <- colSums(dlmest[,i,])
    c(pred.at[i], "mean" = mean(cs), quantile(cs, ci.lims))
  })))
  colnames(cumexp) <- c("vals", "mean", "lower", "upper")

  # Matrix of DLNM surface means and CIs, and plot data
  plot.dat <- as.data.frame(matrix(0, (Lags * length(pred.at)), 10))
  matfit <- matrix(0, length(pred.at), Lags)
  cilower <- matrix(0, length(pred.at), Lags)
  ciupper <- matrix(0, length(pred.at), Lags)
  rownames(matfit) <- rownames(cilower) <- rownames(ciupper) <- pred.at
  colnames(plot.dat) <-  c("Tmin", "Tmax", "Xmin", "Xmax", "PredVal", "Est", "SD", "CIMin", "CIMax", "Effect")
  for (i in 1:Lags) {
    for (j in 1:length(pred.at)) {
      if (object$monotone & zirt.cond) {
        if (sum(splitIter[,i]) == 0 | splitProb[i] < 0.5) {
          plot.dat[(i - 1) * length(pred.at) + j, ] <-
            c(i - 1, i, edge.vals[j], edge.vals[j + 1], pred.at[j],
              0, 0, 0, 0, 0)
          next;
        }
        coordest <- dlmest[i,j,which(splitIter[,i] == 1)]
      } else {
        coordest <- dlmest[i,j,]  
      }
      me <- mean(coordest)
      s <- sd(coordest)
      ci <- quantile(coordest, ci.lims)
      effect <- ifelse(min(ci) > 0, 1, ifelse(max(ci) < 0, -1, 0))
      plot.dat[(i - 1) * length(pred.at) + j, ] <-
        c(i - 1, i, edge.vals[j], edge.vals[j + 1], pred.at[j],
          me, s, ci, effect)
      matfit[j, i] <- mean(coordest)
      cilower[j, i] <- ci[1]
      ciupper[j, i] <- ci[2]
    }
  }
  plot.dat$Effect <- factor(plot.dat$Effect, levels = c(-1, 0, 1), labels = c("-", " ", "+"))

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
              "plot.dat" = plot.dat,
              "matfit" = matfit,
              "cilower" = cilower,
              "ciupper" = ciupper,
              "cenval" = cenval,
              "cumulative.effect" = cumexp,
              "pred.at" = pred.at,
              "gamma.mean" = gamma.mean,
              "gamma.ci" = gamma.ci,
              "logBF" = logBF,
              "splitProb" = splitProb,
              "splitIter" = splitIter)
  class(ret) <- "summary.tdlnm"
  return(ret)
}
