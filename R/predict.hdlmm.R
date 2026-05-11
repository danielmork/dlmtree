#' @method predict hdlmm
#' @rdname predict
predict.hdlmm <- function(x, 
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

  if (!all(x$modNames %in% colnames(new.data))) {
    stop("`new.data` must have the same colunm names as the original model")
  }

  if (!all(sapply(new.exposure.data, function(mat) {is.matrix(mat) && all(is.numeric(mat))}))) {
    stop("`new.exposure.data` must be a list of numeric matrices")
  }

  if (!all(sapply(new.exposure.data, function(mat) {nrow(mat) == nrow(new.data)}))) {
    stop("`new.data` and `new.exposure.data` matrices must have same number of rows")
  }

  if (!all(sapply(new.exposure.data, function(mat) {ncol(mat) == x$pExp}))) {
    stop("`new.data` and `new.exposure.data` matrices must have same number of rows")
  }
    
  ci.lims <- c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2)

  # ---- Setup covariate matrix and modifiers ----
  new.mf  <- model.frame(delete.response(x$formula), 
                          new.data,
                          na.action = na.fail, 
                          drop.unused.levels = TRUE)#,
                          #xlev = x$termLevels)
  z       <- model.matrix(delete.response(x$formula), 
                          new.mf)#,
                          #xlev = x$termLevels)

  n   <- nrow(z)
  mod <- lapply(x$modNames, function(m) {
    if (!is.numeric(x$MoUnique[[m]])) {
      if (!all(unique(new.data[[m]]) %in% x$MoUnique[[m]])) {
        stop("column ", m, " of `new.data` contains additional categories not in original dataset")
      }
    } else {
      if (min(new.data[[m]]) < x$MoUnique[[m]][1] || max(new.data[[m]]) > x$MoUnique[[m]][2]) {
        if (verbose) {
          warning("column ", m, "of `new.data` exceeds range of original dataset")
        }
      }
    }
    new.data[[m]]
    }
  )
  names(mod) <- x$modNames

  # ---- Predict z^T * gamma ----
  ztg.draws     <- z %*% t(x$gamma)
  out$ztg       <- rowMeans(ztg.draws)
  out$ztg.lims  <- apply(ztg.draws, 1, quantile, probs = ci.lims)

  # ---- Predict DLMs ----
  # Main effect
  mark_main <- ceiling(nrow(x$TreeStructs) * x$nExp / 42)
  if (verbose) {
    cat(paste0("Reanalyzing trees for new.data (main): % complete\n",
              "[0--------25--------50--------75--------100]\n '"))
  }
    
  # Main effect for mixture
  main_draws  <- list()
  draws       <- lapply(1:x$mcmcIter, function(i) matrix(0.0, n, x$pExp)) 
  for (exp in x$expNames) {
    main_draws[[exp]] <- draws
  }

  # Iterate through MCMC to add the estimate
  if (is.null(x$fixedIdx)) { 
    for (exp in x$expNames) { 
      for (i in 1:nrow(x$TreeStructs)) {
        # progress bar
        if (verbose && (i %% mark_main == 0)) {
          cat("'")
        }

        # Check if it is the selected exposure
        if ((x$TreeStructs$exp[i] + 1) != which(x$expNames == exp)) {
          next
        }

        Iter <- x$TreeStructs$Iter[i]
        rule <- x$TreeStructs$Rule[i]
        if (rule == "") {
          idx <- 1:n
        } else {
          idx <- which(eval(parse(text = rule)))
        }

        t   <- x$TreeStructs$tmin[i]:x$TreeStructs$tmax[i]
        est <- x$TreeStructs$est[i]
        main_draws[[exp]][[Iter]][idx, t] <- main_draws[[exp]][[Iter]][idx, t] + est
      }
    }
  } else { # Specified MCMC iterations
    for (exp in x$expNames) { # Extra step here for mixture setting
      for (i in 1:nrow(x$TreeStructs)) {
        if (verbose && i %% mark_main == 0) {
          cat("'")
        }

        Iter  <- x$TreeStructs$Iter[i]
        idx   <- fixed.idx[[x$TreeStructs$fixedIdx[i] + 1]]

        t     <- x$TreeStructs$tmin[i]:x$TreeStructs$tmax[i]
        est   <- x$TreeStructs$est[i]
        main_draws[[exp]][[Iter]][idx, t] <- main_draws[[exp]][[Iter]][idx, t] + est
      }
    }
  }

  main_draws <- lapply(main_draws, function(exp) {array(do.call(c, exp), c(n, x$pExp, x$mcmcIter))})

  # Raw predicted DLM and interval
  if (est.dlm) {
    if (verbose) {
      cat("\nEstimating individualized DLM main effects...")
    }
    
    # Output data structure
    out$dlmest <- vector("list", length = x$nExp)
    for (exp in x$expNames) {
      out$dlmest[[exp]] <- list()
      out$dlmest[[exp]][["dlmest"]]       <- sapply(1:x$pExp, function(t) {rowMeans(main_draws[[exp]][,t,,drop=FALSE])}) 
      out$dlmest[[exp]][["dlmest.lower"]] <- sapply(1:x$pExp, function(t) {apply(main_draws[[exp]][,t,,drop=FALSE], 1, quantile, probs = 0.025)})
      out$dlmest[[exp]][["dlmest.upper"]] <- sapply(1:x$pExp, function(t) {apply(main_draws[[exp]][,t,,drop=FALSE], 1, quantile, probs = 0.975)})
    }
  }

  if (x$interaction != 0) {
    # Interaction effect for a mixture setting
    mark_mix <- ceiling(nrow(x$MIX) / 42)
    if (verbose) {
      cat(paste0("\nReanalyzing trees for new.data (interaction): % complete\n",
                "[0--------25--------50--------75--------100]\n '"))
    }

    # Exp i x Exp j (i > j) (mix_mcmc)
    # mix_draws[[Exp x Exp]][lag x lag x n x mcmc]
    mix_draws <- list()
    draws     <- lapply(1:x$mcmcIter, function(i) array(0.0, dim = c(x$pExp, x$pExp, n))) 
    for (mix in x$mixNames) {
      mix_draws[[mix]] <- draws
    }

    mixCount = 1
    for (i in sort(unique(x$MIX$exp1))) {
      for (j in sort(unique(x$MIX$exp2))) {
        # Subsetting matrix with exposure combination
        mixname <- x$mixNames[mixCount]
        newMIX  <- x$MIX[which(x$MIX$exp1 == i & x$MIX$exp2 == j), ]

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

    # Convert from a list form of [[MCMC]](p x p x n) to (p x p x n x MCMC) 4D array form
    mix_draws <- lapply(mix_draws, function(draws) { 
      array(do.call(c, draws), c(x$pExp, x$pExp, n, x$mcmcIter))
    })

    # Return raw estimate of interaction effect
    if (est.dlm) {
      if (verbose) {
        cat("\nEstimating individualized DLM interaction effects...")
      }

      # 4D array calculation components
      matMean <- function(array) {apply(array, c(1, 2), mean)} # 3D array matrix slice mean
      matQt   <- function(mat, qt) {apply(mat, c(1, 2), quantile, probs = qt)}
      grid    <- expand.grid(1:x$pExp, 1:x$pExp)
      
      # Output data structure
      out$mixest <- vector("list", length = length(x$mixNames))
      for (mix in x$mixNames) {
        out$mixest[[mix]] <- vector("list", length = n)
        for (i in 1:n) {
          out$mixest[[mix]][[i]]$mixest       <- matrix(mapply(function(x, y) {matMean(mix_draws[[mix]][x,y,i,,drop=FALSE])}, c(grid$Var1), c(grid$Var2)), nrow = x$pExp)
          out$mixest[[mix]][[i]]$mixest.lower <- matrix(mapply(function(x, y) {matQt(mix_draws[[mix]][x,y,i,,drop=FALSE], 0.025)}, c(grid$Var1), c(grid$Var2)), nrow = x$pExp)
          out$mixest[[mix]][[i]]$mixest.upper <- matrix(mapply(function(x, y) {matQt(mix_draws[[mix]][x,y,i,,drop=FALSE], 0.975)}, c(grid$Var1), c(grid$Var2)), nrow = x$pExp)
        }
      }
    }
  }


  if (verbose) {
    cat("\nCalculating predicted response...")
  }
    
  # Estimate fhat
  fhat.draws <- matrix(0, nrow = n, ncol = x$mcmcIter)

  # Main effect
  for (exp in x$expNames) {
    fhat.draws <- fhat.draws + do.call(rbind, lapply(1:n, function(i) {t(t(main_draws[[exp]][i,,]) %*% new.exposure.data[[exp]][i,]) }))
  }

  if (x$interaction != 0) {
    # Interaction effect
    for (mix in names(mix_draws)) {
      e1  <- unlist(strsplit(mix, "-"))[[1]] # First exposure name
      e2  <- unlist(strsplit(mix, "-"))[[2]] # Second exposure name

      tmp <- do.call(rbind, lapply(1:n, function(i) { # (1 x p)(p x p)(p x 1) per MCMC
            unlist(lapply(lapply(1:x$mcmcIter, function(j) {t(new.exposure.data[[e2]][i, ]) %*% mix_draws[[mix]][,,i,j]}), function(iter) {iter %*% new.exposure.data[[e1]][i, ]}))
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
    y.draws     <- ztg.draws + fhat.draws + sapply(1:x$mcmcIter, function(i) rnorm(nrow(ztg.draws), 0, sqrt(x$sigma2[i])))
    out$y       <- rowMeans(y.draws)
    out$y.lims  <- apply(y.draws, 1, quantile, probs = ci.lims)

    return(out)
  } else if (type == "waic") {
    if (length(outcome) != nrow(ztg.draws)) {
      stop("must provide complete data outcome to calculate WAIC")
    }

    probs <- sapply(1:x$mcmcIter, function(i) dnorm(outcome, ztg.draws[,i] + fhat.draws[,i], sqrt(x$sigma2[i])))
    LPD   <- sum(log(rowMeans(probs)))
    pwaic <- sum(apply(log(probs), 1, var))
    return(list(waic = -2 * (LPD - pwaic), 
                LPD = LPD, 
                pwaic = pwaic)
          )
  }
}
