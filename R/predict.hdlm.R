#' @method predict hdlm
#' @rdname predict
predict.hdlm <- function(x, 
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
    
  new.exposure.data <- matrix(new.exposure.data, ncol = x$pExp)

  if (!is.data.frame(new.data)) {
    stop("`new.data` must be a data.frame with same colunm names as original model")
  }

  if (!all(x$modNames %in% colnames(new.data))) {
    stop("`new.data` must have the same colunm names as the original model")
  }
    
  if (!is.numeric(new.exposure.data) || !is.matrix(new.exposure.data)) {
    stop("`new.exposure.data` must be a numeric matrix")
  }
    
  if (nrow(new.data) != nrow(new.exposure.data)) {
    stop("`new.data` and `new.exposure.data` must have same number of rows")
  }
    
  if (ncol(new.exposure.data) != x$pExp) {
    stop("`new.exposure.data` must have the same number of measurments (columns) as original exposure data")
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
  mark <- ceiling(nrow(x$TreeStructs) / 42)
  if (verbose) {
    cat(paste0("\nReanalyzing trees for new.data: % complete\n",
               "[0--------25--------50--------75--------100]\n '"))
  }

  draws <- lapply(1:x$mcmcIter, function(i) matrix(0.0, n, x$pExp)) # Creates a list of which elements are empty n x p matrix

  # Iterate through MCMC to add the estimate
  if (is.null(x$fixedIdx)) { # No specified MCMC iterations
    for (i in 1:nrow(x$TreeStructs)) {
      if (verbose && (i %% mark == 0)) {
        cat("'")
      }
      Iter <- x$TreeStructs$Iter[i] # Iteration
      rule <- x$TreeStructs$Rule[i] # Rule
      if (rule == "") { 
        idx <- 1:n
      } else { 
        idx <- which(eval(parse(text = rule)))
      }
    }
  } else { # When we specify MCMC iterations
    for (i in 1:nrow(x$TreeStructs)) {
      if (verbose && i %% mark == 0) {
        cat("'")
      }
       
      Iter  <- x$TreeStructs$Iter[i]
      idx   <- fixed.idx[[x$TreeStructs$fixedIdx[i] + 1]]

    }
  }

  # Generate mcmcIter number of matrices (n x pExp)
  draws <- array(do.call(c, draws), c(n, x$pExp, x$mcmcIter)) 

  if (est.dlm) {
    if (verbose) {
      cat("\nestimating individualized DLMs...")
    }
      
    out$dlmest <- sapply(1:x$pExp, function(t) { 
      rowMeans(draws[,t,,drop=FALSE])
    })
    out$dlmest.lower <- sapply(1:x$pExp, function(t) {
      apply(draws[,t,,drop=FALSE], 1, quantile, probs = 0.025)
    })
    out$dlmest.upper <- sapply(1:x$pExp, function(t) {
      apply(draws[,t,,drop=FALSE], 1, quantile, probs = 0.975)
    })
  }
  
  if (verbose) {
    cat("\nestimating exposure effects...")
  }


    fhat.draws <- do.call(rbind, lapply(1:n, function(i) {
    t(t(draws[i,,]) %*% new.exposure.data[i,])
  }))
  out$fhat      <- rowMeans(fhat.draws) # fhat mean for all observation
  out$fhat.lims <- apply(fhat.draws, 1, quantile, probs = ci.lims) # fhat quantiles for all observation



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