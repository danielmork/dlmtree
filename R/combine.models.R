#' combine.models
#'
#' @description Method for combining information from DLMs of single exposure
#'
#' @param mlist a list of models
#'
#' @return A data frame with model fit information of the models included in the list
#' @export combine.models
#'
combine.models <- function(mlist) {
  out     <- mlist[[1]]
  iter    <- out$mcmcIter
  out$DLM <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$DLM
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  out$mcmcIter      <- mlist[[1]]$mcmcIter * length(mlist)
  out$nIter         <- mlist[[1]]$nIter * length(mlist)
  colnames(out$DLM) <- colnames(mlist[[1]]$DLM)
  out$sigma2        <- do.call(c, lapply(mlist, function(l) l$sigma2))
  out$kappa         <- do.call(c, lapply(mlist, function(l) l$kappa))
  out$nu            <- do.call(c, lapply(mlist, function(l) l$nu))
  out$tau           <- do.call(rbind, lapply(mlist, function(l) l$tau))
  out$termNodes     <- do.call(rbind, lapply(mlist, function(l) l$termNodes))
  out$gamma         <- do.call(rbind, lapply(mlist, function(l) l$gamma))
  if(out$monotone) {
    out$zirtGamma   <- do.call(rbind, lapply(mlist, function(l) l$zirtGamma ))
    out$zirtCov     <- do.call(rbind, lapply(mlist, function(l) l$zirtCov))
    out$timeProbs   <- do.call(rbind, lapply(mlist, function(l) l$timeProbs))
    out$timeCounts  <- do.call(rbind, lapply(mlist, function(l) l$timeCounts))
  }
  out$sigma2        <- do.call(c, lapply(mlist, function(l) l$sigma2))
  # out$Yhat <- rowMeans(do.call(cbind, lapply(mlist, function(l) l$Yhat)))
  return(out)
}


#' combine.models.tdlmm
#'
#' @description Method for combining information from DLMs of mixture exposures.
#'
#' @param mlist a list of models
#'
#' @return A data frame with model fit information of the models included in the list
#' @export combine.models.tdlmm
#'
combine.models.tdlmm <- function(mlist) {
  out           <- mlist[[1]]
  iter          <- out$mcmcIter
  out$mcmcIter  <- mlist[[1]]$mcmcIter * length(mlist)
  out$nIter     <- mlist[[1]]$nIter * length(mlist)

  out$DLM <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$DLM
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  out$MIX <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$MIX
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  colnames(out$DLM) <- colnames(mlist[[1]]$DLM)
  colnames(out$MIX) <- colnames(mlist[[1]]$MIX)

  out$expCount    <- do.call(rbind, lapply(mlist, function(l) l$expCount))
  out$expInf      <- do.call(rbind, lapply(mlist, function(l) l$expInf))
  out$expProb     <- do.call(rbind, lapply(mlist, function(l) l$expProb))
  out$muExp       <- do.call(rbind, lapply(mlist, function(l) l$muExp))

  out$mixCount    <- do.call(rbind, lapply(mlist, function(l) l$mixCount))
  out$mixInf      <- do.call(rbind, lapply(mlist, function(l) l$mixInf))
  out$muMix       <- do.call(rbind, lapply(mlist, function(l) l$muMix))

  out$termNodes   <- do.call(rbind, lapply(mlist, function(l) l$termNodes))
  out$termNodes2  <- do.call(rbind, lapply(mlist, function(l) l$termNodes2))
  out$tree1Exp    <- do.call(rbind, lapply(mlist, function(l) l$tree1Exp))
  out$tree2Exp    <- do.call(rbind, lapply(mlist, function(l) l$tree2Exp))


  out$sigma2      <- do.call(c, lapply(mlist, function(l) l$sigma2))
  out$kappa       <- do.call(c, lapply(mlist, function(l) l$kappa))
  out$nu          <- do.call(c, lapply(mlist, function(l) l$nu))
  out$tau         <- do.call(rbind, lapply(mlist, function(l) l$tau))
  out$gamma       <- do.call(rbind, lapply(mlist, function(l) l$gamma))
  return(out)
}
