#' summary.tdlmm
#'
#' @title Creates a summary object of class 'tdlmm'
#' @description Method for creating a summary object of class 'tdlmm'
#'
#' @param object an object of type 'tdlmm', the output from tdlmm()
#' @param conf.level confidence level (default = 0.95)
#' @param marginalize value(s) for calculating marginal DLMs, defaults to
#' "mean", can also specify a percentile from 1-99 for all other exposures, or
#' a named vector with specific values for each exposure
#' @param log10BF.crit Bayes Factor criteria for selecting exposures and
#' interactions, such that log10(BayesFactor) > x. Default = 0.5
#' @param verbose show progress in console
#' @param keep.mcmc keep all mcmc iterations (large memory requirement)
#' @param ... additional parameters
#'
#' @returns list of type 'summary.tdlmm'
#' @export
#'
summary.tdlmm <- function(object,
                          conf.level = 0.95,
                          marginalize = "mean",
                          log10BF.crit = 0.5,
                          verbose = TRUE,
                          keep.mcmc = FALSE,
                          ...)
{
  res               <- list()
  res$nIter         <- object$nIter
  res$nThin         <- object$nThin
  res$mcmcIter      <- floor(res$nIter / res$nThin)
  res$nBurn         <- object$nBurn
  res$nTrees        <- object$nTrees
  res$treePrior     <- object$treePriorTDLM
  res$nLags         <- max(object$TreeStructs$tmax)
  res$nExp          <- object$nExp
  res$nMix          <- object$nMix
  res$expNames      <- object$expNames
  res$mixNames      <- object$mixNames
  res$interaction   <- object$interaction
  res$mixPrior      <- object$mixPrior
  res$formula       <- object$formula
  res$formula.zi    <- object$formula.zi
  res$family        <- object$family
  res$droppedCovar  <- object$droppedCovar
  res$conf.level    <- conf.level
  res$ci.lims       <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  res$log10BF.crit  <- log10BF.crit


  # ---- Set levels for marginaliztion ----
  if (marginalize[1] == "mean") {
    res$marg.values <- sapply(object$X, function(i) i$intX)
  } else if (length(marginalize) == 1 & is.numeric(marginalize)) {
    if (marginalize < 0 | marginalize > 100) {
      stop("A percentile for marginialization must be between 0 and 100")
    }
    res$marg.values <- sapply(object$X, function(i) {i$Xquant[round(marginalize) + 1] })
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
    BF <- (log10(colMeans(object$expCount > 0)) -
             log10(colMeans(object$expCount == 0))) -
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
    idx <- which(object$TreeStructs$exp == (i - 1))
    est <- dlmEst(as.matrix(object$TreeStructs[idx, -(3:4)]), res$nLags, res$mcmcIter)
    # est: dim 1 = time, dim 2 = iteration
    res$DLM[[i]] <- list("mcmc" = est,
                         "marg" = array(0, dim(est)),
                         "name" = object$expNames[i])
  }
  names(res$DLM)  <- res$expNames
  iqr_plus_mean   <- function(i) c(quantile(i, 0.25), mean(i), quantile(i, 0.75))
  res$expInc      <- apply(object$expCount > 0, 2, mean)

  if (res$nExp == 1) {
    res$expVar <- apply(object$mixCount > 0, 2, mean)
    #names(res$expVar) <- names(res$expInc)
  } else {
    res$expVar <- apply((apply(object$muExp * object$expInf, 1, rank) - 1),
                                1, iqr_plus_mean) / (res$nExp - 1)
  }


  # ---- Mixture MCMC samples ----
  if (verbose) {
    cat("Reconstructing interaction effects...\n0%...")
  }
    
  res$MIX <- list()
  nMix    <- 1
  for (i in sort(unique(object$MIX$exp1))) {
    for (j in sort(unique(object$MIX$exp2))) {
      if (verbose) {
        if (nMix == ceiling(0.25 * res$nMix)) {cat("25%...")}
        if (nMix == ceiling(0.5 * res$nMix)) {cat("50%...")}
        if (nMix == ceiling(0.75 * res$nMix)) {cat("75%...")}
        if (nMix == res$nMix) {cat("100%\n")}
        nMix <- nMix + 1
      }

      idx <- which(object$MIX$exp1 == i & object$MIX$exp2 == j)
      if (length(idx) > 0) {
        est <- mixEst(as.matrix(object$MIX[idx,,drop = FALSE]), res$nLags, res$mcmcIter)
        m   <- paste0(object$expNames[i + 1], "-", object$expNames[j + 1])
        res$MIX[[m]] <-
          list("matfit"   =  sapply(1:res$nLags, function(k) rowMeans(est[,k,,drop=FALSE])),
               "cilower"  =  sapply(1:res$nLags, function(k) {
                              apply(est[,k,], 1, quantile, probs = res$ci.lims[1]) }),
               "ciupper"  =  sapply(1:res$nLags, function(k) {
                              apply(est[,k,], 1, quantile, probs = res$ci.lims[2]) }),
               "rows"     = object$expNames[i + 1],
               "cols"     = object$expNames[j + 1])

        if (keep.mcmc) {
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
          
          if (keep.mcmc) {
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
    res$mixInc <- apply(object$mixCount > 0, 2, mean)
    if (res$nMix == 1) {
      res$mixVar <- c(1)
      names(res$mixVar) <- names(res$mixInc)
    } else {
      res$mixVar <- apply((apply(object$muMix * object$mixInf, 1, rank) - 1),
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
    if (keep.mcmc) {
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

    # DLM non-linear effects
    # if (any(names(res$MIX) == paste0(ex.name, "-", ex.name))) {
    #   quad <- sapply(1:res$mcmcIter, function(i) {
    #     diag(res$MIX[[paste0(ex.name, "-", ex.name)]]$mcmc[,,i]) })
    #   med <- object$X[[ex.name]]$Xquant["50%"] *
    #     (marg - object$X[[ex.name]]$intX * quad +
    #        object$X[[ex.name]]$Xquant["50%"] * quad)
    #   res$DLM[[ex.name]]$marg.nonlin.matfit <-
    #     sapply(object$X[[ex.name]]$Xquant, function(x) {
    #       apply(x * (marg - object$X[[ex.name]]$intX * quad + x * quad) - med,
    #             1, mean)
    #     })
    #   res$DLM[[ex.name]]$marg.nonlin.cilower <-
    #     sapply(object$X[[ex.name]]$Xquant, function(x) {
    #       apply(x * (marg - object$X[[ex.name]]$intX * quad + x * quad) - med,
    #             1, quantile, probs = res$ci.lims[1])
    #     })
    #   res$DLM[[ex.name]]$marg.nonlin.ciupper <-
    #     sapply(object$X[[ex.name]]$Xquant, function(x) {
    #       apply(x * (marg - object$X[[ex.name]]$intX * quad + x * quad) - med,
    #             1, quantile, probs = res$ci.lims[2])
    #     })
    #   res$DLM[[ex.name]]$marg.nonlin.cw <-
    #     (res$DLM[[ex.name]]$marg.nonlin.cilower > 0 |
    #        res$DLM[[ex.name]]$marg.nonlin.ciupper < 0)
    #
    #   colnames(res$DLM[[ex.name]]$marg.nonlin.matfit) <-
    #     colnames(res$DLM[[ex.name]]$marg.nonlin.cilower) <-
    #     colnames(res$DLM[[ex.name]]$marg.nonlin.ciupper) <-
    #     colnames(res$DLM[[ex.name]]$marg.nonlin.cw) <-
    #     as.character(object$X[[ex.name]]$Xquant)
    # }
  }

  # ---- Fixed effect estimates ----
  if (verbose) {
    cat("Calculating fixed effects...\n")
  }
  if (!object$zinb) {
    # Gaussian / Logistic
    res$gamma.mean  <- colMeans(object$gamma)
    res$gamma.ci    <- apply(object$gamma, 2, quantile, probs = res$ci.lims)
  } else {
    # ZINB
    # binary
    res$b1.mean <- colMeans(object$b1)
    res$b1.ci   <- apply(object$b1, 2, quantile, probs = res$ci.lims)

    # count
    res$b2.mean <- colMeans(object$b2)
    res$b2.ci   <- apply(object$b2, 2, quantile, probs = res$ci.lims)

    # Dispersion parameter
    res$r.mean  <- mean(object$r)
    res$r.ci    <- quantile(object$r, probs = res$ci.lims)
  }

  # ---- Return ----
  res$sig2noise <- ifelse(is.null(object$sigma2), NA,
                          var(object$fhat) / mean(object$sigma2))
  res$rse       <- sd(object$sigma2)
  res$n         <- nrow(object$data)

  class(res) <- "summary.tdlmm"
  
  return(res)
}
