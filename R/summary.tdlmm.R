#' @method summary tdlmm
#' @rdname summary
#'
#' @export
summary.tdlmm <- function(x, conf.level = 0.95, marginalize = "mean", log10BF.crit = 0.5, mcmc = FALSE, verbose = TRUE, ...)
{
  ret <- list(
    "ctr" = list(class    = x$class,
                 n.trees  = x$nTrees,
                 n.iter   = x$nIter,
                 n.thin   = x$nThin,
                 n.burn   = x$nBurn,
                 response = x$family),
    "conf.level"    = conf.level,
    "n.lag"         = max(x$TreeStructs$tmax),
    "n.exp"         = x$nExp,
    "n.mix"         = x$nMix,
    "exp.names"     = x$expNames,
    "mix.names"     = x$mixNames,
    "interaction"   = x$interaction,
    "mix.prior"     = x$mixPrior,
    "formula"       = x$formula,
    "dropped.covar" = x$droppedCovar,
    "ci.lims"       = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2),
    "log10BF.crit"  = log10BF.crit
  )
  
  mcmcIter     <- floor(x$nIter / x$nThin)

  # ---- Set levels for marginalization ----
  if (marginalize[1] == "mean") {
    ret$marg.values <- sapply(x$X, function(i) i$intX)
  } else if (length(marginalize) == 1 & is.numeric(marginalize)) {
    if (marginalize < 0 | marginalize > 100) {
      stop("A percentile for marginialization must be between 0 and 100")
    }
    ret$marg.values <- sapply(x$X, function(i) {i$Xquant[round(marginalize) + 1] })
  } else if (length(marginalize) == ret$n.exp & is.numeric(marginalize)) {
    ret$marg.values <- marginalize
  } else {
    stop("`marginalize` is incorrectly specified, see ?summary.tdlmm for details")
  }
  names(ret$marg.values) <- ret$exp.names


  if(verbose){
    cat("Specified co-exposure values:\n")
    for(exp in 1:length(ret$exp.names)){
      cat("-", ret$exp.names[exp], ":", ret$marg.values[exp])
      cat("\n")
    }
    cat("\n")
  }

  # ---- Bayes factor variable selection ----
  ret$expSel        <- rep(FALSE, ret$n.exp)
  names(ret$expSel) <- ret$exp.names
  if (ret$mix.prior > 0) {
    # individual exposures
    priorProbInc <- 1 - exp(lgamma(2*ret$ctr$n.trees + ret$n.exp*ret$mix.prior) +
                          lgamma((ret$n.exp + 1)*ret$mix.prior) -
                          lgamma(2*ret$ctr$n.trees + (ret$n.exp + 1)*ret$mix.prior) -
                          lgamma(ret$n.exp*ret$mix.prior))
    BF <- (log10(colMeans(x$expCount > 0)) -
             log10(colMeans(x$expCount == 0))) -
      (log10(priorProbInc) - log10(1 - priorProbInc))

    ret$expBF   <- 10^BF
    ret$expSel  <- (BF > log10BF.crit)
  }

  # ---- Main effect MCMC samples ----
  if (verbose) {
    cat("Reconstructing main effects...\n")
  }
  ret$DLM <- list()
  for (i in 1:ret$n.exp) {
    idx <- which(x$TreeStructs$exp == (i - 1))
    est <- dlmEst(as.matrix(x$TreeStructs[idx, -(3:4)]), ret$n.lag, mcmcIter)
    # est: dim 1 = time, dim 2 = iteration
    ret$DLM[[i]] <- list("mcmc" = est,
                         "marg" = array(0, dim(est)),
                         "name" = x$expNames[i])
  }
  names(ret$DLM)  <- ret$exp.names
  iqr_plus_mean   <- function(i) c(quantile(i, 0.25), mean(i), quantile(i, 0.75))
  ret$expInc      <- apply(x$expCount > 0, 2, mean)

  if (ret$n.exp == 1) {
    ret$expVar <- apply(x$mixCount > 0, 2, mean)
    #names(ret$expVar) <- names(ret$expInc)
  } else {
    ret$expVar <- apply((apply(x$muExp * x$expInf, 1, rank) - 1),
                                1, iqr_plus_mean) / (ret$n.exp - 1)
  }


  # ---- Mixture MCMC samples ----
  if (verbose) {
    cat("Reconstructing interaction effects...\n0%...")
  }
    
  ret$MIX <- list()
  nMix    <- 1
  for (i in sort(unique(x$MIX$exp1))) {
    for (j in sort(unique(x$MIX$exp2))) {
      if (verbose) {
        if (nMix == ceiling(0.25 * ret$n.mix)) {cat("25%...")}
        if (nMix == ceiling(0.5 * ret$n.mix)) {cat("50%...")}
        if (nMix == ceiling(0.75 * ret$n.mix)) {cat("75%...")}
        if (nMix == ret$n.mix) {cat("100%\n")}
        nMix <- nMix + 1
      }

      idx <- which(x$MIX$exp1 == i & x$MIX$exp2 == j)
      if (length(idx) > 0) {
        est <- mixEst(as.matrix(x$MIX[idx,,drop = FALSE]), ret$n.lag, mcmcIter)
        m   <- paste0(x$expNames[i + 1], "-", x$expNames[j + 1])
        ret$MIX[[m]] <-
          list("matfit"   =  sapply(1:ret$n.lag, function(k) rowMeans(est[,k,,drop=FALSE])),
               "cilower"  =  sapply(1:ret$n.lag, function(k) {
                              apply(est[,k,], 1, quantile, probs = ret$ci.lims[1]) }),
               "ciupper"  =  sapply(1:ret$n.lag, function(k) {
                              apply(est[,k,], 1, quantile, probs = ret$ci.lims[2]) }),
               "rows"     = x$expNames[i + 1],
               "cols"     = x$expNames[j + 1])

        if (mcmc) {
          ret$MIX[[m]]$mcmc <- est
        }

        # Calculate marginal effects
        ret$DLM[[i + 1]]$marg <- ret$DLM[[i + 1]]$marg +
          t(sapply(1:ret$n.lag, function(k) colSums(est[k,,]))) *
          ret$marg.values[j + 1] * ifelse(i == j, 0.5, 1)

        ret$DLM[[j + 1]]$marg <- ret$DLM[[j + 1]]$marg +
          t(sapply(1:ret$n.lag, function(k) colSums(est[,k,]))) *
          ret$marg.values[i + 1] * ifelse(i == j, 0.5, 1)


        # Fold surface of self interaction
        if (i == j) {
          est <- 0.5 * est * array(upper.tri(diag(ret$n.lag), diag = TRUE), dim(est)) +
            0.5 * aperm(est, c(2, 1, 3)) *
            array(upper.tri(diag(ret$n.lag), diag = TRUE), dim(est))
          
          if (mcmc) {
            ret$MIX[[m]]$mcmc   <- est
          }

          ret$MIX[[m]]$matfit   <- sapply(1:ret$n.lag, function(k) {rowMeans(est[,k,,drop=FALSE]) })
          ret$MIX[[m]]$cilower  <- sapply(1:ret$n.lag, function(k) {apply(est[,k,], 1,  quantile, probs = ret$ci.lims[1]) })
          ret$MIX[[m]]$ciupper  <- sapply(1:ret$n.lag, function(k) {apply(est[,k,], 1, quantile, probs = ret$ci.lims[2]) })
        }

        ret$MIX[[m]]$cw <- (ret$MIX[[m]]$cilower > 0 | ret$MIX[[m]]$ciupper < 0)

        # Range of confidence levels for plots
        mixCIs <- lapply(1:ret$n.lag, function(k) {apply(est[,k,], 1, quantile, probs = c(1:10/200, 190:199/200))})

        ciProbs <- c(99:90/100, 0)
        ret$MIX[[m]]$cw.plot <-
          sapply(1:ret$n.lag, function(k) {
            ciProbs[
              sapply(1:ret$n.lag, function(l) {
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

  if (ret$interaction > 0) {
    ret$mixInc <- apply(x$mixCount > 0, 2, mean)
    if (ret$n.mix == 1) {
      ret$mixVar <- c(1)
      names(ret$mixVar) <- names(ret$mixInc)
    } else {
      ret$mixVar <- apply((apply(x$muMix * x$mixInf, 1, rank) - 1),
                          1, iqr_plus_mean) / (ret$n.mix - 1)
    }
  }

  # ---- DLM marginal effects ----
  if (verbose) {
    cat("Calculating marginal effects...\n")
  }

  for (ex.name in names(ret$DLM)) {

    # DLM marginal effects
    marg <- ret$DLM[[ex.name]]$mcmc + ret$DLM[[ex.name]]$marg
    if (mcmc) {
      ret$DLM[[ex.name]]$marg <- marg
    } else {
      ret$DLM[[ex.name]]$mcmc <- NULL
      ret$DLM[[ex.name]]$marg <- NULL
    }

    ret$DLM[[ex.name]]$marg.matfit  <- apply(marg, 1, mean)
    ret$DLM[[ex.name]]$marg.cilower <- apply(marg, 1, quantile, probs = ret$ci.lims[1])
    ret$DLM[[ex.name]]$marg.ciupper <- apply(marg, 1, quantile, probs = ret$ci.lims[2])
    ret$DLM[[ex.name]]$marg.cw      <- (ret$DLM[[ex.name]]$marg.cilower > 0 | ret$DLM[[ex.name]]$marg.ciupper < 0)

    # Cumulative effects
    cumulative <- colSums(marg)
    ret$DLM[[ex.name]]$cumulative <-
      list("mean"     = mean(cumulative),
           "ci.lower" = quantile(cumulative, ret$ci.lims[1]),
           "ci.upper" = quantile(cumulative, ret$ci.lims[2]))

  }

  # ---- Fixed effect estimates ----
  if (verbose) {
    cat("Calculating fixed effects...\n")
  }
  
  if (!x$zinb) {
    # Gaussian / Logistic
    ret$gamma.mean  <- colMeans(x$gamma)
    ret$gamma.ci    <- apply(x$gamma, 2, quantile, probs = ret$ci.lims)
  } else {
    # ZINB
    ret$formula.zi  <- x$formula.zi
    
    # binary
    ret$b1.mean <- colMeans(x$b1)
    ret$b1.ci   <- apply(x$b1, 2, quantile, probs = ret$ci.lims)

    # count
    ret$b2.mean <- colMeans(x$b2)
    ret$b2.ci   <- apply(x$b2, 2, quantile, probs = ret$ci.lims)

    # Dispersion parameter
    ret$r.mean  <- mean(x$r)
    ret$r.ci    <- quantile(x$r, probs = ret$ci.lims)
  }

  # ---- Return ----
  ret$sig.to.noise <- ifelse(is.null(x$sigma2), NA,
                          var(x$fhat) / mean(x$sigma2))
  ret$rse       <- sd(x$sigma2)
  ret$n         <- nrow(x$data)

  class(ret) <- "summary.tdlmm"
  
  return(ret)
}
