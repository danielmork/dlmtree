combine.models <- function(mlist) {
  out <- mlist[[1]]
  iter <- out$mcmcIter
  out$DLM <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$DLM
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  out$mcmcIter <- mlist[[1]]$mcmcIter * length(mlist)
  out$nIter <- mlist[[1]]$nIter * length(mlist)
  colnames(out$DLM) <- colnames(mlist[[1]]$DLM)
  out$sigma2 <- do.call(c, lapply(mlist, function(l) l$sigma2))
  out$kappa <- do.call(c, lapply(mlist, function(l) l$kappa))
  out$nu <- do.call(c, lapply(mlist, function(l) l$nu))
  out$tau <- do.call(rbind, lapply(mlist, function(l) l$tau))
  out$termNodes <- do.call(rbind, lapply(mlist, function(l) l$termNodes))
  out$gamma <- do.call(rbind, lapply(mlist, function(l) l$gamma))
  out$zirt <- do.call(rbind, lapply(mlist, function(l) l$zirt))
  out$zirtCov <- do.call(rbind, lapply(mlist, function(l) l$zirtCov))
  out$timeProbs <- do.call(rbind, lapply(mlist, function(l) l$timeProbs))
  out$timeCounts <- do.call(rbind, lapply(mlist, function(l) l$timeCounts))
  out$sigma2 <- do.call(c, lapply(mlist, function(l) l$sigma2))
  out$Yhat <- rowMeans(do.call(cbind, lapply(mlist, function(l) l$Yhat)))
  return(out)
}
