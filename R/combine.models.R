#' combine.models
#' 
#' @title Combines information from DLMs of single exposure
#' @description Method for combining information from DLMs of single exposure
#'
#' @param mlist a list of models
#'
#' @returns A data frame with model fit information of the models included in the list
#' @export combine.models
#'
combine.models <- function(mlist) {
  out     <- mlist[[1]]
  if (length(unique(sapply(mlist, function(l) l$mcmcIter))) != 1)
    stop("all chains must have same number of posterior samples")
  iter    <- out$mcmcIter
  out$TreeStructs <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$TreeStructs
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  colnames(out$TreeStructs) <- colnames(mlist[[1]]$TreeStructs)
  out$mcmcIter      <- mlist[[1]]$mcmcIter * length(mlist)
  out$nIter         <- mlist[[1]]$nIter * length(mlist)
  out$sigma2        <- do.call(c, lapply(mlist, function(l) l$sigma2))
  out$kappa         <- do.call(c, lapply(mlist, function(l) l$kappa))
  out$nu            <- do.call(c, lapply(mlist, function(l) l$nu))
  out$tau           <- do.call(rbind, lapply(mlist, function(l) l$tau))
  out$termNodes     <- do.call(rbind, lapply(mlist, function(l) l$termNodes))
  out$gamma         <- do.call(rbind, lapply(mlist, function(l) l$gamma))
  out$fhat          <- rowMeans(do.call(cbind, lapply(mlist, function(l) l$fhat)))
  out$Yhat          <- rowMeans(do.call(cbind, lapply(mlist, function(l) l$Yhat)))
  out$sigma2        <- do.call(c, lapply(mlist, function(l) l$sigma2))
  if(out$monotone) {
    out$zirtGamma   <- do.call(rbind, lapply(mlist, function(l) l$zirtGamma ))
    out$timeProbs   <- do.call(rbind, lapply(mlist, function(l) l$timeProbs))
    out$zirtSplitCounts  <- do.call(rbind, lapply(mlist, function(l) l$zirtSplitCounts))
  }
  return(out)
}


#' combine.models.tdlmm
#'
#' @title Combines information from DLMs of mixture exposures.
#' @description Method for combining information from DLMs of mixture exposures.
#'
#' @param mlist a list of models
#'
#' @returns A data frame with model fit information of the models included in the list
#' @export combine.models.tdlmm
#'
combine.models.tdlmm <- function(mlist) {
  out           <- mlist[[1]]
  iter          <- out$mcmcIter
  out$mcmcIter  <- mlist[[1]]$mcmcIter * length(mlist)
  out$nIter     <- mlist[[1]]$nIter * length(mlist)

  out$TreeStructs <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$TreeStructs
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  out$MIX <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$MIX
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  colnames(out$TreeStructs) <- colnames(mlist[[1]]$TreeStructs)
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
