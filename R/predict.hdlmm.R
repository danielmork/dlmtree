#' predict.hdlmm
#'
#' @title Calculates predicted response for HDLMM
#' @description Method for calculating predicted response for HDLMM
#'
#' @param object fitted dlmtree model with class hdlmm
#' @param new.data new data frame which contains the same covariates and modifiers used to fit HDLMM model
#' @param new.exposure.data new data frame/list which contains the same length of exposure lags used to fit HDLMM model
#' @param ci.level credible interval level for posterior predictive distribution
#' @param type type of prediction: "response" (default) or "waic". "waic" must be specified with `outcome` parameter
#' @param outcome outcome required for WAIC calculation
#' @param fixed.idx fixed index
#' @param est.dlm flag for estimating dlm effect
#' @param verbose TRUE (default) or FALSE: print output
#' @param ... additional parameters
#'
#' @returns Posterior predictive distribution draws
#' @export
predict.hdlmm <- function(object, 
                         new.data, 
                         new.exposure.data,
                         ci.level = 0.95, 
                         type = "response", 
                         outcome = NULL,
                         fixed.idx = list(), 
                         est.dlm = FALSE, 
                         verbose = TRUE,
                         ...)
{ 
  `%notin%` <- Negate(`%in%`)
  out       <- list()

  if (!is.data.frame(new.data)) {
    stop("`new.data` must be a data.frame with same colunm names as original model")
  }

  if (!all(object$modNames %in% colnames(new.data))) {
    stop("`new.data` must have the same colunm names as the original model")
  }

  if (!all(sapply(new.exposure.data, function(mat) {is.matrix(mat) && all(is.numeric(mat))}))) {
    stop("`new.exposure.data` must be a list of numeric matrices")
  }

  if (!all(sapply(new.exposure.data, function(mat) {nrow(mat) == nrow(new.data)}))) {
    stop("`new.data` and `new.exposure.data` matrices must have same number of rows")
  }

  if (!all(sapply(new.exposure.data, function(mat) {ncol(mat) == object$pExp}))) {
    stop("`new.data` and `new.exposure.data` matrices must have same number of rows")
  }
    
  ci.lims <- c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2)

  # ---- Setup covariate matrix and modifiers ----
  new.mf  <- model.frame(delete.response(object$formula), 
                          new.data,
                          na.action = na.fail, 
                          drop.unused.levels = TRUE,
                          xlev = object$termLevels)
  z       <- model.matrix(delete.response(object$formula), 
                          new.mf,
                          xlev = object$termLevels)

  n   <- nrow(z)
  mod <- lapply(object$modNames, function(m) {
    if (!is.numeric(object$MoUnique[[m]])) {
      if (!all(unique(new.data[[m]]) %in% object$MoUnique[[m]])) {
        stop("column ", m, " of `new.data` contains additional categories not in original dataset")
      }
    } else {
      if (min(new.data[[m]]) < object$MoUnique[[m]][1] || max(new.data[[m]]) > object$MoUnique[[m]][2]) {
        if (verbose) {
          warning("column ", m, "of `new.data` exceeds range of original dataset")
        }
      }
    }
    new.data[[m]]
    }
  )
  names(mod) <- object$modNames

  # ---- Predict z^T * gamma ----
  ztg.draws     <- z %*% t(object$gamma)
  out$ztg       <- rowMeans(ztg.draws)
  out$ztg.lims  <- apply(ztg.draws, 1, quantile, probs = ci.lims)

  # ---- Predict DLMs ----
  # Main effect
  mark_main <- ceiling(nrow(object$TreeStructs) * object$nExp / 42)
  if (verbose) {
    cat(paste0("Reanalyzing trees for new.data (main): % complete\n",
              "[0--------25--------50--------75--------100]\n '"))
  }
    
  # Main effect for mixture
  main_draws  <- list()
  draws       <- lapply(1:object$mcmcIter, function(i) matrix(0.0, n, object$pExp)) # Creates a list of which elements are empty n x p matrix    
  for (exp in object$expNames) {
    main_draws[[exp]] <- draws
  }

  # Iterate through MCMC to add the estimate
  if (is.null(object$fixedIdx)) { # No specified MCMC iterations
    for (exp in object$expNames) { # Extra step here for mixture setting
      for (i in 1:nrow(object$TreeStructs)) {
        # progress bar
        if (verbose && (i %% mark_main == 0)) {
          cat("'")
        }

        # Check if it is the selected exposure
        if ((object$TreeStructs$exp[i] + 1) != which(object$expNames == exp)) {
          next
        }

        Iter <- object$TreeStructs$Iter[i]
        rule <- object$TreeStructs$Rule[i]
        if (rule == "") {
          idx <- 1:n
        } else {
          idx <- which(eval(parse(text = rule)))
        }

        t   <- object$TreeStructs$tmin[i]:object$TreeStructs$tmax[i]
        est <- object$TreeStructs$est[i]
        main_draws[[exp]][[Iter]][idx, t] <- main_draws[[exp]][[Iter]][idx, t] + est
      }
    }
  } else { # Specified MCMC iterations
    for (exp in object$expNames) { # Extra step here for mixture setting
      for (i in 1:nrow(object$TreeStructs)) {
        if (verbose && i %% mark_main == 0) {
          cat("'")
        }

        Iter  <- object$TreeStructs$Iter[i]
        idx   <- fixed.idx[[object$TreeStructs$fixedIdx[i] + 1]]

        t     <- object$TreeStructs$tmin[i]:object$TreeStructs$tmax[i]
        est   <- object$TreeStructs$est[i]
        main_draws[[exp]][[Iter]][idx, t] <- main_draws[[exp]][[Iter]][idx, t] + est
      }
    }
  }

  # do.call generates c(draws)
  # Convert from a list form of [[MCMC]](nxp) to (n x p x MCMC) array form
  main_draws <- lapply(main_draws, function(exp) {array(do.call(c, exp), c(n, object$pExp, object$mcmcIter))})

  # Raw predicted DLM and interval
  if (est.dlm) {
    if (verbose) {
      cat("\nEstimating individualized DLM main effects...")
    }
    
    # Output data structure
    out$dlmest <- vector("list", length = object$nExp)
    for (exp in object$expNames) {
      out$dlmest[[exp]] <- list()
      out$dlmest[[exp]][["dlmest"]]       <- sapply(1:object$pExp, function(t) {rowMeans(main_draws[[exp]][,t,,drop=FALSE])}) 
      out$dlmest[[exp]][["dlmest.lower"]] <- sapply(1:object$pExp, function(t) {apply(main_draws[[exp]][,t,,drop=FALSE], 1, quantile, probs = 0.025)})
      out$dlmest[[exp]][["dlmest.upper"]] <- sapply(1:object$pExp, function(t) {apply(main_draws[[exp]][,t,,drop=FALSE], 1, quantile, probs = 0.975)})
    }
  }

  if (object$interaction != 0) {
    # Interaction effect for a mixture setting
    mark_mix <- ceiling(nrow(object$MIX) / 42)
    if (verbose) {
      cat(paste0("\nReanalyzing trees for new.data (interaction): % complete\n",
                "[0--------25--------50--------75--------100]\n '"))
    }

    # Exp i x Exp j (i > j) (mix_mcmc)
    # mix_draws[[Exp x Exp]][lag x lag x n x mcmc]
    mix_draws <- list()
    draws     <- lapply(1:object$mcmcIter, function(i) array(0.0, dim = c(object$pExp, object$pExp, n))) 
    for (mix in object$mixNames) {
      mix_draws[[mix]] <- draws
    }

    mixCount = 1
    for (i in sort(unique(object$MIX$exp1))) {
      for (j in sort(unique(object$MIX$exp2))) {
        # Subsetting matrix with exposure combination
        mixname <- object$mixNames[mixCount]
        newMIX  <- object$MIX[which(object$MIX$exp1 == i & object$MIX$exp2 == j), ]

        if (i == j) { # HDLMM only consider no-self interaction
          next
        }
        
        # MCMC draws with rules
        for (k in 1:nrow(newMIX)) {
          if (verbose && k %% mark_mix == 0)
            cat("'")

          Iter <- newMIX$Iter[k]
          Rule <- newMIX$Rule[k]

          if (Rule == "") {
            idx <- 1:n
          } else {
            idx <- which(eval(parse(text = Rule)))
            if (length(idx) == 0) {
              idx = 0
            } # avoid integer(0) problem when n = 1 for individualized dlm
          }
          
          t1range <- newMIX$tmin1[k]:newMIX$tmax1[k]
          t2range <- newMIX$tmin2[k]:newMIX$tmax2[k]
          est     <- newMIX$est[k]

          for (t1 in t1range) {
            for (t2 in t2range) {
              mix_draws[[mixname]][[Iter]][t1, t2, idx] <- mix_draws[[mixname]][[Iter]][t1, t2, idx] + est
            }
          }
        }

        mixCount = mixCount + 1
      }
    }

    # do.call generates c(draws)
    # Convert from a list form of [[MCMC]](p x p x n) to (p x p x n x MCMC) 4D array form
    # for each pairwise interaction
    mix_draws <- lapply(mix_draws, function(draws) { 
      array(do.call(c, draws), c(object$pExp, object$pExp, n, object$mcmcIter))
    })

    # Return raw estimate of interaction effect
    if (est.dlm) {
      if (verbose) {
        cat("\nEstimating individualized DLM interaction effects...")
      }

      # 4D array calculation components
      matMean <- function(array) {apply(array, c(1, 2), mean)} # 3D array matrix slice mean
      matQt   <- function(mat, qt) {apply(mat, c(1, 2), quantile, probs = qt)}
      grid    <- expand.grid(1:object$pExp, 1:object$pExp)
      
      # Output data structure
      out$mixest <- vector("list", length = length(object$mixNames))
      for (mix in object$mixNames) {
        out$mixest[[mix]] <- vector("list", length = n)
        for (i in 1:n) {
          out$mixest[[mix]][[i]]$mixest       <- matrix(mapply(function(x, y) {matMean(mix_draws[[mix]][x,y,i,,drop=FALSE])}, c(grid$Var1), c(grid$Var2)), nrow = object$pExp)
          out$mixest[[mix]][[i]]$mixest.lower <- matrix(mapply(function(x, y) {matQt(mix_draws[[mix]][x,y,i,,drop=FALSE], 0.025)}, c(grid$Var1), c(grid$Var2)), nrow = object$pExp)
          out$mixest[[mix]][[i]]$mixest.upper <- matrix(mapply(function(x, y) {matQt(mix_draws[[mix]][x,y,i,,drop=FALSE], 0.975)}, c(grid$Var1), c(grid$Var2)), nrow = object$pExp)
        }
      }
    }
  }


  if (verbose) {
    cat("\nCalculating predicted response...")
  }
    
  # draws[i,,]: extract all MCMC samples of ith observation where col of a matrix is the MCMC sample: p x MCMC
  # t(draws[i,,]): MCMC x p
  # (t(draws[i,,]) %*% new.exposure.data[i,]): (MCMC x p)x(p x 1)
  # t(t(draws[i,,]) %*% new.exposure.data[i,]): 1 x MCMC
  # do.call(rbind, lapply) -> Do the above for all observations and combine them to a data frame using rbind -> n x MCMC
  # Estimate fhat
  fhat.draws <- matrix(0, nrow = n, ncol = object$mcmcIter)

  # Main effect
  for (exp in object$expNames) {
    fhat.draws <- fhat.draws + do.call(rbind, lapply(1:n, function(i) {t(t(main_draws[[exp]][i,,]) %*% new.exposure.data[[exp]][i,]) }))
  }

  if (object$interaction != 0) {
    # Interaction effect
    for (mix in names(mix_draws)) {
      e1  <- unlist(strsplit(mix, "-"))[[1]] # First exposure name
      e2  <- unlist(strsplit(mix, "-"))[[2]] # Second exposure name

      tmp <- do.call(rbind, lapply(1:n, function(i) { # (1 x p)(p x p)(p x 1) per MCMC
            unlist(lapply(lapply(1:object$mcmcIter, function(j) {t(new.exposure.data[[e2]][i, ]) %*% mix_draws[[mix]][,,i,j]}), function(iter) {iter %*% new.exposure.data[[e1]][i, ]}))
            # left = lapply(1:object$mcmcIter, function(j) {t(new.exposure.data[[e2]][i, ]) %*% mix_draws[[mix]][,,i,j]})
            # right = lapply(left, function(iter) {iter %*% new.exposure.data[[e1]][i, ]})
            # return(unlist(right))
          }
        )
      )

      fhat.draws <- fhat.draws + tmp
    }

    rm(main_draws, mix_draws)
  } else {
    rm(main_draws)
  }
  

  # Final result
  out$fhat.draws  <- fhat.draws
  out$fhat        <- rowMeans(fhat.draws)
  out$fhat.lims   <- apply(fhat.draws, 1, quantile, probs = ci.lims)  
  

  # ---- Outcome predictions ----
  # fixed effect + DLM + normal error
  if (type == "response") {
    y.draws     <- ztg.draws + fhat.draws + sapply(1:object$mcmcIter, function(i) rnorm(nrow(ztg.draws), 0, sqrt(object$sigma2[i])))
    out$y       <- rowMeans(y.draws)
    out$y.lims  <- apply(y.draws, 1, quantile, probs = ci.lims)

    return(out)
  } else if (type == "waic") {
    if (length(outcome) != nrow(ztg.draws)) {
      stop("must provide complete data outcome to calculate WAIC")
    }

    probs <- sapply(1:object$mcmcIter, function(i) dnorm(outcome, ztg.draws[,i] + fhat.draws[,i], sqrt(object$sigma2[i])))
    LPD   <- sum(log(rowMeans(probs)))
    pwaic <- sum(apply(log(probs), 1, var))
    return(list(waic = -2 * (LPD - pwaic), 
                LPD = LPD, 
                pwaic = pwaic)
          )
  }
}
