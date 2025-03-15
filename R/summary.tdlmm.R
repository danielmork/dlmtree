#' @method summary tdlmm
#' @rdname summary
#'
#' @export
summary.tdlmm <- function(x, conf.level = 0.95, marginalize = "mean", log10BF.crit = 0.5, mcmc = FALSE, verbose = TRUE, ...)
{
  res               <- list()
  res$nIter         <- x$nIter
  res$nThin         <- x$nThin
  res$mcmcIter      <- floor(res$nIter / res$nThin)
  res$nBurn         <- x$nBurn
  res$nTrees        <- x$nTrees
  res$treePrior     <- x$treePriorTDLM
  res$nLags         <- max(x$TreeStructs$tmax)
  res$nExp          <- x$nExp
  res$nMix          <- x$nMix
  res$expNames      <- x$expNames
  res$mixNames      <- x$mixNames
  res$interaction   <- x$interaction
  res$mixPrior      <- x$mixPrior
  res$formula       <- x$formula
  res$formula.zi    <- x$formula.zi
  res$family        <- x$family
  res$droppedCovar  <- x$droppedCovar
  res$conf.level    <- conf.level
  res$ci.lims       <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  res$log10BF.crit  <- log10BF.crit


  # ---- Set levels for marginaliztion ----
  if (marginalize[1] == "mean") {
    res$marg.values <- sapply(x$X, function(i) i$intX)
  } else if (length(marginalize) == 1 & is.numeric(marginalize)) {
    if (marginalize < 0 | marginalize > 100) {
      stop("A percentile for marginialization must be between 0 and 100")
    }
    res$marg.values <- sapply(x$X, function(i) {i$Xquant[round(marginalize) + 1] })
  } else if (length(marginalize) == res$nExp & is.numeric(marginalize)) {
    res$marg.values <- marginalize
  } else {
    stop("`marginalize` is incorrectly specified, see ?summary.tdlmm for details")
  }
  names(res$marg.values) <- res$expNames


  if(verbose){
    cat("Specified co-exposure values:\n")
    for(exp in 1:length(res$expNames)){
      cat("-", res$expNames[exp], ":", res$marg.values[exp])
      cat("\n")
    }
    cat("\n")
  }

  # ---- Bayes factor variable selection ----
  res$expSel        <- rep(FALSE, res$nExp)
  names(res$expSel) <- res$expNames
  if (res$mixPrior > 0) {
    # individual exposures
    priorProbInc <- 1 - exp(lgamma(2*res$nTrees + res$nExp*res$mixPrior) +
                          lgamma((res$nExp + 1)*res$mixPrior) -
                          lgamma(2*res$nTrees + (res$nExp + 1)*res$mixPrior) -
                          lgamma(res$nExp*res$mixPrior))
    BF <- (log10(colMeans(x$expCount > 0)) -
             log10(colMeans(x$expCount == 0))) -
      (log10(priorProbInc) - log10(1 - priorProbInc))

    res$expBF   <- 10^BF
    res$expSel  <- (BF > log10BF.crit)
  }

  # ---- Main effect MCMC samples ----
  if (verbose) {
    cat("Reconstructing main effects...\n")
  }
  res$DLM <- list()
  for (i in 1:res$nExp) {
    idx <- which(x$TreeStructs$exp == (i - 1))
    est <- dlmEst(as.matrix(x$TreeStructs[idx, -(3:4)]), res$nLags, res$mcmcIter)
    # est: dim 1 = time, dim 2 = iteration
    res$DLM[[i]] <- list("mcmc" = est,
                         "marg" = array(0, dim(est)),
                         "name" = x$expNames[i])
  }
  names(res$DLM)  <- res$expNames
  iqr_plus_mean   <- function(i) c(quantile(i, 0.25), mean(i), quantile(i, 0.75))
  res$expInc      <- apply(x$expCount > 0, 2, mean)

  if (res$nExp == 1) {
    res$expVar <- apply(x$mixCount > 0, 2, mean)
    #names(res$expVar) <- names(res$expInc)
  } else {
    res$expVar <- apply((apply(x$muExp * x$expInf, 1, rank) - 1),
                                1, iqr_plus_mean) / (res$nExp - 1)
  }


  # ---- Mixture MCMC samples ----
  if (verbose) {
    cat("Reconstructing interaction effects...\n0%...")
  }
    
  res$MIX <- list()
  nMix    <- 1
  for (i in sort(unique(x$MIX$exp1))) {
    for (j in sort(unique(x$MIX$exp2))) {
      if (verbose) {
        if (nMix == ceiling(0.25 * res$nMix)) {cat("25%...")}
        if (nMix == ceiling(0.5 * res$nMix)) {cat("50%...")}
        if (nMix == ceiling(0.75 * res$nMix)) {cat("75%...")}
        if (nMix == res$nMix) {cat("100%\n")}
        nMix <- nMix + 1
      }

      idx <- which(x$MIX$exp1 == i & x$MIX$exp2 == j)
      if (length(idx) > 0) {
        est <- mixEst(as.matrix(x$MIX[idx,,drop = FALSE]), res$nLags, res$mcmcIter)
        m   <- paste0(x$expNames[i + 1], "-", x$expNames[j + 1])
        res$MIX[[m]] <-
          list("matfit"   =  sapply(1:res$nLags, function(k) rowMeans(est[,k,,drop=FALSE])),
               "cilower"  =  sapply(1:res$nLags, function(k) {
                              apply(est[,k,], 1, quantile, probs = res$ci.lims[1]) }),
               "ciupper"  =  sapply(1:res$nLags, function(k) {
                              apply(est[,k,], 1, quantile, probs = res$ci.lims[2]) }),
               "rows"     = x$expNames[i + 1],
               "cols"     = x$expNames[j + 1])

        if (mcmc) {
          res$MIX[[m]]$mcmc <- est
        }

        # Calculate marginal effects
        res$DLM[[i + 1]]$marg <- res$DLM[[i + 1]]$marg +
          t(sapply(1:res$nLags, function(k) colSums(est[k,,]))) *
          res$marg.values[j + 1] * ifelse(i == j, 0.5, 1)

        res$DLM[[j + 1]]$marg <- res$DLM[[j + 1]]$marg +
          t(sapply(1:res$nLags, function(k) colSums(est[,k,]))) *
          res$marg.values[i + 1] * ifelse(i == j, 0.5, 1)


        # Fold surface of self interaction
        if (i == j) {
          est <- 0.5 * est * array(upper.tri(diag(res$nLags), diag = TRUE), dim(est)) +
            0.5 * aperm(est, c(2, 1, 3)) *
            array(upper.tri(diag(res$nLags), diag = TRUE), dim(est))
          
          if (mcmc) {
            res$MIX[[m]]$mcmc   <- est
          }

          res$MIX[[m]]$matfit   <- sapply(1:res$nLags, function(k) {rowMeans(est[,k,,drop=FALSE]) })
          res$MIX[[m]]$cilower  <- sapply(1:res$nLags, function(k) {apply(est[,k,], 1,  quantile, probs = res$ci.lims[1]) })
          res$MIX[[m]]$ciupper  <- sapply(1:res$nLags, function(k) {apply(est[,k,], 1, quantile, probs = res$ci.lims[2]) })
        }

        res$MIX[[m]]$cw <- (res$MIX[[m]]$cilower > 0 | res$MIX[[m]]$ciupper < 0)

        # Range of confidence levels for plots
        mixCIs <- lapply(1:res$nLags, function(k) {apply(est[,k,], 1, quantile, probs = c(1:10/200, 190:199/200))})

        ciProbs <- c(99:90/100, 0)
        res$MIX[[m]]$cw.plot <-
          sapply(1:res$nLags, function(k) {
            ciProbs[
              sapply(1:res$nLags, function(l) {
                min(c(11,which(
                  sapply(1:10, function(p) {
                    (mixCIs[[k]][p, l] > 0 | mixCIs[[k]][21 - p, l] < 0)
                  })
                )))
              })
            ]
          })
      }
    }
  }

  if (res$interaction > 0) {
    res$mixInc <- apply(x$mixCount > 0, 2, mean)
    if (res$nMix == 1) {
      res$mixVar <- c(1)
      names(res$mixVar) <- names(res$mixInc)
    } else {
      res$mixVar <- apply((apply(x$muMix * x$mixInf, 1, rank) - 1),
                          1, iqr_plus_mean) / (res$nMix - 1)
    }
  }

  # ---- DLM marginal effects ----
  if (verbose) {
    cat("Calculating marginal effects...\n")
  }

  for (ex.name in names(res$DLM)) {

    # DLM marginal effects
    marg <- res$DLM[[ex.name]]$mcmc + res$DLM[[ex.name]]$marg
    if (mcmc) {
      res$DLM[[ex.name]]$marg <- marg
    } else {
      res$DLM[[ex.name]]$mcmc <- NULL
      res$DLM[[ex.name]]$marg <- NULL
    }

    res$DLM[[ex.name]]$marg.matfit  <- apply(marg, 1, mean)
    res$DLM[[ex.name]]$marg.cilower <- apply(marg, 1, quantile, probs = res$ci.lims[1])
    res$DLM[[ex.name]]$marg.ciupper <- apply(marg, 1, quantile, probs = res$ci.lims[2])
    res$DLM[[ex.name]]$marg.cw      <- (res$DLM[[ex.name]]$marg.cilower > 0 | res$DLM[[ex.name]]$marg.ciupper < 0)

    # Cumulative effects
    cumulative <- colSums(marg)
    res$DLM[[ex.name]]$cumulative <-
      list("mean"     = mean(cumulative),
           "ci.lower" = quantile(cumulative, res$ci.lims[1]),
           "ci.upper" = quantile(cumulative, res$ci.lims[2]))

  }

  # ---- Fixed effect estimates ----
  if (verbose) {
    cat("Calculating fixed effects...\n")
  }
  if (!x$zinb) {
    # Gaussian / Logistic
    res$gamma.mean  <- colMeans(x$gamma)
    res$gamma.ci    <- apply(x$gamma, 2, quantile, probs = res$ci.lims)
  } else {
    # ZINB
    # binary
    res$b1.mean <- colMeans(x$b1)
    res$b1.ci   <- apply(x$b1, 2, quantile, probs = res$ci.lims)

    # count
    res$b2.mean <- colMeans(x$b2)
    res$b2.ci   <- apply(x$b2, 2, quantile, probs = res$ci.lims)

    # Dispersion parameter
    res$r.mean  <- mean(x$r)
    res$r.ci    <- quantile(x$r, probs = res$ci.lims)
  }

  # ---- Return ----
  res$sig2noise <- ifelse(is.null(x$sigma2), NA,
                          var(x$fhat) / mean(x$sigma2))
  res$rse       <- sd(x$sigma2)
  res$n         <- nrow(x$data)

  class(res) <- "summary.tdlmm"
  
  return(res)
}
