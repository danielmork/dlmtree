#' summary.tdlnm
#'
#' @param object an object of class 'tdlnm', result of a call to tdlnm
#' @param pred.at numerical vector of exposure values to make predictions for
#' at each time period
#' @param cenval scalar exposure value that acts as a reference point for
#' predictions at all other exposure values
#' @param conf.level confidence level for computation of credible intervals
#' @param exposure.se smoothing factor, if different from model
#'
#' @return
#' @export
#'
#' @examples
summary.tdlnm <- function(object, 
                          pred.at = c(), 
                          cenval = 0, 
                          conf.level = 0.95,
                          exposure.se = NULL)
{
  Iter <-  max(object$DLM$Iter)
  Trees <- max(object$DLM$Tree)
  Lags <- length(object$Mo$Tsplits) + 1
  Xvals <- c(object$Mo$X)
  n <- length(object$Y)
  ci.lims <- c((1-conf.level)/2, 1-(1-conf.level)/2)
  if (is.null(exposure.se) & !is.null(object$Mo$SE))
    exposure.se <- mean(object$Mo$SE)
  else if (is.null(exposure.se))
    exposure.se = 0

  # Determine how to break up data for summary
  if (is.null(pred.at)) {
    Xsplits <- object$Mo$Xsplits
    if (length(Xsplits) == 0)
      Xsplits <- seq(min(Xvals), max(Xvals), length.out = 102)
    edge.quants <- sort(unique(c(0, ecdf(Xvals)(Xsplits), 1)))
    Xsplits <- quantile(Xvals, edge.quants)
    half.quants <- rowMeans(cbind(edge.quants[-1],edge.quants[-length(edge.quants)]))
    pred.vals <- quantile(Xvals, half.quants)
  } else {
    half.quants <- sort(unique(ecdf(Xvals)(pred.at)))
    pred.vals <- quantile(Xvals, half.quants)
    edge.quants <- unique(c(0, rowMeans(cbind(half.quants[-1],half.quants[-length(half.quants)])), 1))
    Xsplits <- quantile(Xvals, edge.quants)
  }

  # Define actual viewing locations as middle quantile between splits
  cen.val <- cenval
  cat("Centered DLNM at exposure value", cenval, "\n")
  if (exposure.se == 0) {
    cen.quant <- which.min(abs(pred.vals - cenval))
    if (length(object$Mo$Xsplits) == 0)
      dlmest <- dlnmEst(as.matrix(object$DLM), Xsplits, Lags, Iter, cen.quant, exposure.se, FALSE, TRUE)
    else
      dlmest <- dlnmEst(as.matrix(object$DLM), Xsplits, Lags, Iter, cen.quant, exposure.se, FALSE, FALSE)
  } else {
    dlmest <- dlnmEst(as.matrix(object$DLM), Xsplits, Lags, Iter, cen.val, exposure.se, TRUE, FALSE)
  }

  # DLNM Estimates
  cumexp <- as.data.frame(t(sapply(1:length(half.quants), function(i) {
    cs <- colSums(dlmest[,i,])
    c("mean" = mean(cs), quantile(cs, ci.lims))
  })))
  dlmmean <- as.data.frame(matrix(0, (Lags * (length(half.quants))), 9))
  matfit <- matrix(0, length(half.quants), Lags)
  cilower <- matrix(0, length(half.quants), Lags)
  ciupper <- matrix(0, length(half.quants), Lags)
  rownames(matfit) <- rownames(cilower) <- rownames(ciupper) <- pred.vals
  colnames(dlmmean) <-  c("Tmin", "Tmax", "Xmin", "Xmax", "Est", "SD", "CIMin", "CIMax", "Effect")
  for (i in 1:Lags) {
    for (j in 1:(length(half.quants))) {
      coordest <- dlmest[i,j,]
      m <- mean(coordest)
      s <- sd(coordest)
      ci <- quantile(coordest, ci.lims)
      effect <- ifelse(min(ci) > 0, 1, ifelse(max(ci) < 0, -1, 0))
      dlmmean[(i - 1) * (length(Xsplits) - 1) + j, ] <-
        c(i - 1, i, Xsplits[j], Xsplits[j + 1],
          m, s, ci, effect)
      matfit[j, i] <- mean(coordest)
      cilower[j, i] <- ci[1]
      ciupper[j, i] <- ci[2]
    }
  }
  dlmmean$Effect <- factor(dlmmean$Effect, levels = c(-1, 0, 1), labels = c("-", " ", "+"))
  
  # Fixed effect estimates
  gamma.mean <- colMeans(object$gamma)
  gamma.ci <- apply(object$gamma, 2, quantile, probs = ci.lims)

  # Return
  ret <- list("ctr" = object$ctr,
              "conf.level" = conf.level,
              "sig.to.noise" = var(object$fhat) / mean(object$sigma2),
              "dlnm.estimates" = dlmmean, 
              "matfit" = matfit,
              "cilower" = cilower, 
              "ciupper" = ciupper,
              "Xsplits" = Xsplits, 
              "cenval" = cen.val,
              "cumexp" = cumexp, 
              "pred.vals" = pred.vals,
              "halfQuants" = half.quants,
              "gamma.mean" = gamma.mean,
              "gamma.ci" = gamma.ci)
  class(ret) <- "summary.tdlnm"
  return(ret)
}