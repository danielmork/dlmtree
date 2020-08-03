#' tdlnm
#' @title Treed Distributed Lag Non-Linear Models
#' @description TDLNM is a method for estimating a Treed Distirubted Lag
#' Non-Linear Model. It operates by building an ensemble of regression trees,
#' which each contribute a partial estimate of the exposure-time surface using
#' the 'exposure.data', while controlling for a model given by 'formula'.
#'
#' @param formula object of class formula, a symbolic description of the fixed
#' effect model to be fitted, e.g. y ~ a + b
#' @param data data frome containing variables used in the formula
#' @param exposure.data numerical matrix of complete exposure data with same
#' length as data
#' @param exposure.splits scalar indicating the number of splits (divided
#' evenly across quantiles of the exposure data) or list with two components:
#' 'type' set to 'values' or 'quantiles', and 'split.vals' set as a numerical
#' vector indicating the corresponding exposure values or quantiles for splits.
#' Setting exposure.splits equal to 0 will change the model to a distributed lag
#' model, which assumes a linear effect of exposure.
#' @param exposure.se numerical matrix of exposure standard errors with same
#' length as data or scalar smoothing factor representing a uniform standard
#' error over each exposure measurement
#' @param n.trees integer for number of trees in ensemble
#' @param n.burn integer for length of burn-in
#' @param n.iter integer for number of iterations to run model after burn-in
#' @param n.thin integer thinning factor, i.e. keep every tenth iteration
#' @param response.type 'gaussian' (default), 'binary', or 'binomial', if
#' 'binomial', must set binomial.size
#' @param binomial.size integer vector of binomial size for 'binomial' response
#' @param mixture.interactions 'noself' (default) which allows interactions to
#' occur only between two different exposures, 'all' which also allows
#' interactions within the same exposure, or 'none' which eliminates all
#' interactions and estimates only main effects of each exposure
#' @param tree.params numerical vector of alpha and beta hyperparameters
#' controlling tree depth (see Bayesian CART, 1998)
#' @param step.prob numerical vector for probability of 1) grow/prune, 2) change
#' @param time.prior.period integer vector corresponding to columns of exposure.data
#' for which there is prior belief of the beginning/end of a critical window;
#' for example, if a critical window is believed to exists between periods
#' 10 and 15, enter c(10, 15)
#' @param time.prior.certainty integer from 1 to 10 corresponding to prior
#' certainty (10 is absolute certainty) that a critical window exists in
#' 'time.prior.period'
#' @param time.split.prior numeric vector of splitting probabilities for time,
#' must be length 1 less than columns of exposure.data; this entry will nullify
#' entries into time.prior.period and time.prior.certainty
#' @param subset integer vector to analyze only a subset of data and exposures
#' @param verbose true/false print output
#' @param diagnostics true/false keep model diagnostic such as terminal nodes,
#' acceptance details, etc.
#' @param ... NA
#'
#' @usage Model is recommended to be run for at minimum 5000 burn-in iterations
#' followed by 15000 sampling iterations with a thinning factor of 10.
#' Convergence can be checked by re-running the model and validating consistency
#' of results.
#'
#' @return obect of class 'tdlnm' or 'tdlm'
#' @export
#'
tdlnm <- function(formula,
                  data,
                  exposure.data,
                  exposure.splits = 50,
                  exposure.se = NULL,
                  n.trees = 20,
                  n.burn = 2000,
                  n.iter = 5000,
                  n.thin = 10,
                  response.type = "gaussian",
                  binomial.size = NULL,
                  mixture.interactions = "noself",
                  tree.params = c(.95, 2),
                  step.prob = c(.4, .2, .2),
                  time.prior.period = NULL,
                  time.prior.certainty = NULL,
                  time.split.prior = NULL,
                  subset = NULL,
                  verbose = TRUE,
                  diagnostics = FALSE,
                  max.threads = 2,
                  ...)
{
  model <- new.env()
  # ---- Check inputs ----
  options(stringsAsFactors = F)

  # data for formula
  if (!is.data.frame(data))
    stop("`data` must be a data.frame")

  # exposure data (single exposure)
  if (!is.list(exposure.data)) {
    mixture <- FALSE
    if (!is.numeric(exposure.data))
      stop("`exposure.data` must be a numeric matrix")
    if (nrow(data) != nrow(exposure.data))
      stop("`data` and `exposure.data` must have same number of observations")
    # exposure smoothing/error
    if (!is.null(exposure.se)) {
      if (!is.numeric(exposure.se))
        stop("`exposure.se` must be a scalar or numeric matrix")
      if (length(exposure.se) == 1)
        exposure.se <- matrix(exposure.se, nrow(exposure.data), ncol(exposure.data))
      if (nrow(exposure.se) != nrow(exposure.data) &
          ncol(exposure.se) != ncol(exposure.data))
        stop("`exposure.se` and `exposure.data` must have same dimensions")
     }

  # exposure data (multiple exposures)
  } else {
    mixture <- TRUE
    for (i in 1:length(exposure.data)) {
      if (!is.numeric(exposure.data[[i]]))
        stop("each exposure in list `exposure.data` must be a numeric matrix")
      if (nrow(data) != nrow(exposure.data[[i]]))
        stop("`data` and `exposure.data` must have same number of observations")
    }
  }

  # iteration control
  if (all(sapply(list(n.trees,n.burn,n.iter,n.thin),
                 function(i) is.integer(i) & i > 0)))
    stop("n.* must be integer and > 0")

  # response type
  if (!(response.type %in% c("gaussian", "binomial", "binary")))
    stop("`response.type` must be one of `gaussian`, `binomial`, or `binary`")

  # binomial size
  if (response.type == "binomial") {
    if (length(binomial.size) == 1)
      binomial.size <- rep(binomial.size, nrow(data))
    if (length(binomial.size) != nrow(data))
      stop("`binomial.size` must be positive integer and same length as data")
    model$binomialSize <- force(binomial.size)
  } else if (response.type == "binary") {
    binomial.size <- rep(1, nrow(data))
    model$binomialSize <- force(binomial.size)
  }

  # mixture interactions
  if (mixture) {
    if (!(mixture.interactions %in% c("noself", "all", "marginal", "none")))
      stop("`mixture.interactions must be one of `noself`, `all`, `marginal` or `none`")
    if (mixture.interactions %in% c("marginal", "none"))
      model$interaction = 0
    else if (mixture.interactions == "noself")
      model$interaction = 1
    else
      model$interaction = 2
  }

  # tree parameters
  if (length(tree.params) != 2) {
    stop("tree.params must have length 2")
  } else {
    if (tree.params[1] > 1 | tree.params[1] < 0)
      stop("tree.params[1] must be between 0 and 1")
    if (tree.params[2] < 0)
      stop("tree.params[2] must be greater than 0")
  }

  # step probabilities
  if (any(step.prob < 0) | any(step.prob > 1))
    stop("step.prob components must be between zero and one")



  # ---- Setup control and response variables ----
  model$formula <- force(as.formula(formula))
  tf <- terms.formula(model$formula, data = data)
  if (!attr(tf, "response"))
    stop("no valid response in formula")
  model$intercept <- force(ifelse(attr(tf, "intercept"), TRUE, FALSE))
  if (length(which(attr(tf, "term.labels") %in% colnames(data))) == 0)
    stop("no valid variables in formula")




  # ---- Model control arguments ----
  model$nTrees <- force(n.trees)
  model$nBurn <- force(n.burn)
  model$nIter <- force(n.iter)
  model$nThin <- force(n.thin)
  model$response <- force(response.type)
  model$verbose <- force(verbose)
  model$diagnostics <- force(diagnostics)
  model$treePrior <- force(tree.params)
  model$threads <- force(max.threads)
  if (mixture)
    model$stepProb <- force(c(step.prob[1], step.prob[1], step.prob[2], step.prob[3]))
  else
    model$stepProb <- force(c(step.prob[1], step.prob[1], step.prob[2]))
  model$stepProb <- force(model$stepProb / sum(model$stepProb))

  if (model$verbose)
    cat("Preparing data...\n")



  # ---- Create data subset ----
  if (!is.null(subset)) {
    if (length(subset) > 1 & is.integer(subset) &
        all(subset > 0) & all(subset <= nrow(data))) {
      data <- data[subset,]
      if (mixture) {
        model$expNames <- names(exposure.data)
        exposure.data <- lapply(exposure.data, function(i) i[subset,])
        names(exposure.data) <- model$expNames
      } else
        exposure.data <- exposure.data[subset,]
      if (!is.null(exposure.se))
        exposure.se <- exposure.se[subset,]
      if (!is.null(binomial.size))
        model$binomialSize <- force(binomial.size[subset])
    } else {
      stop("invalid subset, must be integers within range of data length")
    }
  }



  # ---- Exposure splits ----
  # multiple exposures
  if (mixture) {
    model$dlFunction <- "tdlmm"
    model$splitProb <- as.double(c())
    model$Xsplits <- as.double(c())
    model$nSplits <- 0
    model$X <- list()
    for (i in 1:length(exposure.data)) {
      model$X[[i]] <- force(list(Xscale = sd(exposure.data[[i]]),
                           Xmean = mean(exposure.data[[i]]),
                           X = exposure.data[[i]]))
      model$X[[i]]$X <- force((model$X[[i]]$X - model$X[[i]]$Xmean) / model$X[[i]]$Xscale)
      model$X[[i]]$Tcalc <- force(sapply(1:ncol(model$X[[i]]$X), function(j) {
        rowSums(model$X[[i]]$X[,1:j,drop=F]) }))
      if (any(is.na(model$X)))
        stop("missing values in exposure data")
    }
    model$expNames <- names(exposure.data)
    model$expProb <- force(rep(1/length(model$X), length(model$X)))

  # single exposure
  } else {
    model$X <- force(exposure.data)
    model$Xrange <- force(range(exposure.data))
    if (is.null(exposure.se)) {
      model$smooth <- FALSE
      model$SE <- matrix(0.0, 0, 0)
    } else {
      model$smooth <- TRUE
      model$SE <- force(exposure.se)
    }

    if (length(exposure.splits) == 1) {

      # no splits -- DLM
      if (exposure.splits == 0) {
        model$dlFunction <- "tdlm"
        model$splitProb <- as.double(c())
        model$Xsplits <- as.double(c())
        model$nSplits <- 0
        model$Xscale <- force(sd(model$X))
        model$Xmean <- force(mean(model$X))
        model$X <- force((model$X - model$Xmean) / model$Xscale)
        model$Tcalc <- force(sapply(1:ncol(model$X), function(i) rowSums(model$X[,1:i,drop=F])))

      # splits defined by quantiles of exposure
      } else {
        model$dlFunction <- "tdlnm"
        if (is.list(exposure.splits)) {

        } else {
          model$Xsplits <- force(sort(unique(quantile(model$X, (1:(exposure.splits - 1)) /
                                                        exposure.splits))))
          model$nSplits <- force(length(model$Xsplits))
          model$splitProb <- force(rep(1 / model$nSplits, model$nSplits))
        }
      }

    # splits defined by specific values or quantiles
    } else {
      model$dlFunction <- "tdlnm"
      if (!is.list(exposure.splits))
        stop("`exposure.splits` must be a named list specifying `type` and `split.vals`")
      if (exposure.splits$type == "values") {
        model$Xsplits <- force(sort(unique(exposure.splits$split.vals)))
        model$Xsplits <- force(model$Xsplits[which(model$Xsplits > min(model$X) &
                                               model$Xsplits < max(model$X))])
      } else if (exposure.splits$type == "quantiles") {
        if (any(exposure.splits$split.vals > 1 | exposure.splits$split.vals < 0))
          stop("`exposure.splits$split.vals` must be between zero and one")
        model$Xsplits <- force(sort(unique(quantile(model$X, exposure.splits$split.vals))))
        model$Xsplits <- force(model$Xsplits[which(model$Xsplits > min(model$X) &
                                               model$Xsplits < max(model$X))])
      } else {
        stop("`exposure.splits$type` must be one of `values` or `quantiles`")
      }
      model$nSplits <- force(length(model$Xsplits))
      model$splitProb <- force(rep(1 / model$nSplits, model$nSplits))
    }

    if (length(model$Xsplits) > 0) {
      model$Xscale <- 1#force(sqrt(nrow(model$X)))# * ncol(model$X)
      model$Tcalc <- force(sapply(1:ncol(model$X),
                                  function(i) rep(i / model$Xscale, nrow(model$X))))
      if (model$smooth) {
        model$Xcalc <- force(sapply(model$Xsplits,
                                    function(i) rowSums(pnorm((i - model$X) / model$SE))) /
                               model$Xscale)
      } else {
        model$Xcalc <- force(sapply(model$Xsplits, function(i) rowSums(model$X < i)) /
                               model$Xscale)
      }
    }
  }




  # ---- Setup time prior ----
  if (mixture)
    pX <- ncol(exposure.data[[1]])
  else
    pX <- ncol(exposure.data)
  model$timeProb <- force(rep(1 / (pX - 1), pX - 1))
  # if (!is.list(exposure.data)) {
  #   pX <- ncol(exposure.data)
  # } else {
  #   pX <- ncol(exposure.data[[1]])
  # }
  # if (!is.null(time.split.prior)) {
  #   if (length(time.split.prior) != pX - 1 | any(time.split.prior < 0))
  #     stop(paste0("'time.split.prior' must be >= 0 and of length ", pX - 1))
  #   model$time.prior <- time.split.prior / sum(time.split.prior)
  # } else {
  #   if (is.null(time.prior.period)) {
  #     model$time.prior <- rep(1 / (pX - 1), pX - 1)
  #   } else {
  #     if (length(time.prior.certainty) > 1 | !all(time.prior.certainty %in% 1:10))
  #       stop("'time.prior.certainty' must be an integer from 1 to 10")
  #     time.prior.certainty <- time.prior.certainty * 10
  #     model$time.prior <- rowMeans(
  #       sapply((time.prior.period / pX), function(p) {
  #         dcauchy(1:pX/(pX + 1), p, 1 / (10 * time.prior.certainty))
  #       }))
  #     model$time.prior <- colMeans(rbind(model$time.prior[1:(pX - 1)],
  #                                        model$time.prior[2:pX]))
  #     model$time.prior <- model$time.prior / sum(model$time.prior)
  #   }
  # }



  # ---- Scale data and setup exposures ----
  data <- droplevels(data)
  mf <- model.frame(model$formula, data = data,
                    drop.unused.levels = TRUE,
                    na.action = NULL)
  if (any(is.na(mf)))
    stop("missing values in model data, use `complete.cases()` to subset data")
  model$Y <- force(model.response(mf))
  model$Z <- force(model.matrix(model$formula, data = mf))
  QR <- qr(crossprod(model$Z))
  model$Z <- model$Z[,sort(QR$pivot[seq_len(QR$rank)])]
  model$Z <- force(scaleModelMatrix(model$Z))



  if (model$response == "gaussian") {
    model$Ymean <- sum(range(model$Y))/2
    model$Yscale <- diff(range(model$Y - model$Ymean))
    model$Y <- force(scale(model$Y, center = model$Ymean,
                     scale = model$Yscale))
  } else {
    model$Yscale <- 1
    model$Ymean <- 0
    model$Y <- force(scale(model$Y, center = 0, scale = 1))
  }
  model$Y <- force(c(model$Y))
  model$Zscale <- attr(model$Z, "scaled:scale")
  model$Zmean <- attr(model$Z, "scaled:center")
  model$Znames <- colnames(model$Z)
  model$Z <- force(matrix(model$Z, nrow(model$Z), ncol(model$Z)))
  # if (model$dlFunction == "tdlnm" && !model$intercept) {
  #   stop("tdlnm requires the intercept to be included in the model fixed effect")
  # }



  # ---- Run model ----
  model$try = 0 # rerun if initial draws contain NaN
  model.list <- lapply(ls(envir = model), function(i) model[[i]])
  names(model.list) <- ls(envir = model)
  if (mixture) {
    if (model$response == "gaussian")
      out <- tdlnmMixGaussian(model.list)
    # else
    #   tdlnm.mix.binomial(model)
  } else {
    if (model$response == "gaussian")
      out <- tdlnmGaussian(model.list)
    else
      out <- tdlnmBinomial(model.list)
  }

  for (n in names(out))
    model[[n]] <- out[[n]]

  # ---- Prepare output ----
  model$Y <- model$Y * model$Yscale + model$Ymean
  model$fhat <- model$fhat * model$Yscale
  model$sigma2 <- model$sigma2 * (model$Yscale ^ 2)
  if (model$diagnostics) {
    model$treeAccept <- as.data.frame(model$treeAccept)
    colnames(model$treeAccept) <- c("step", "success", "term", "treeMhr", "mhr")
  }

  # rescale fixed effect estimates
  if (model$intercept) {
    model$gamma[,-1] <- sapply(2:ncol(model$gamma), function(i) {
      model$gamma[,i] * model$Yscale / model$Zscale[i]})
    model$gamma[,1] <- model$gamma[,1] -
      apply(sapply(2:ncol(model$gamma), function(i) {
        model$gamma[,i] * model$Zmean[i]
      }), 1, sum)
  } else {
    model$gamma <- sapply(1:ncol(model$gamma), function(i) {
      model$gamma[,i] * model$Yscale / model$Zscale[i]})
  }
  colnames(model$gamma) <- model$Znames

  # rescale DLM estimates
  if (mixture) {
    model$DLM <- as.data.frame(model$DLM)
    colnames(model$DLM) <- c("Iter", "Tree", "exp", "tmin", "tmax",
                             "est", "kappa")
    model$MIX <- as.data.frame(model$MIX)
    colnames(model$MIX) <- c("Iter", "Tree", "exp1", "tmin1", "tmax1",
                             "exp2", "tmin2", "tmax2", "est", "kappa")
    for (i in 1:length(model$X)) {
      idx <- which(model$DLM$exp == (i - 1))
      if (length(idx) > 0)
        model$DLM$est[idx] <- model$DLM$est[idx] * model$Yscale / model$X[[i]]$Xscale

      for (j in i:length(model$X)) {
        idx <- which(model$MIX$exp1 == (i - 1) & model$MIX$exp2 == (j - 1))
        if (length(idx) > 0)
          model$MIX$est[idx] <- model$MIX$est[idx] * model$Yscale / (model$X[[i]]$Xscale * model$X[[j]]$Xscale)
      }
    }

  } else {
    model$DLM <- as.data.frame(model$DLM)
    colnames(model$DLM) <- c("Iter", "Tree", "xmin", "xmax", "tmin", "tmax", "est", "kappa")
    model$DLM$est <- model$DLM$est * model$Yscale / model$Xscale
    model$DLM$xmin <- sapply(model$DLM$xmin, function(i) {
      if (i == 0) -Inf
      else model$Xsplits[i]
    })
    model$DLM$xmax <- sapply(model$DLM$xmax, function(i) {
      if (i == (length(model$Xsplits) + 1)) Inf
      else model$Xsplits[i]
    })
  }



  # remove data
  # model$X <- NULL
  # model$Z <- NULL

  # Change env to list
  model.out <- lapply(names(model), function(i) model[[i]])
  names(model.out) <- names(model)
  rm(list = ls(envir = model), envir = model)
  rm(model)
  gc()

  if (model.out$dlFunction == "tdlnm")
    class(model.out) <- "tdlnm"
  else if (model.out$dlFunction == "tdlm")
    class(model.out) <- "tdlm"
  else if (model.out$dlFunction == "tdlmm")
    class(model.out) <- "tdlmm"
  return(model.out)
}
