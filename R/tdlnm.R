#' tdlnm
#' @title Treed Distributed Lag Non-Linear Models
#' @description TDLNM is a method for estimating Distributed Lag
#' Linear and Non-Linear Models (DLMs/DLNMs). It operates by building an
#' ensemble of regression trees, which each partition the exposure-time-
#' response surface and make estimates at each sector. Trees from the ensemble
#' each contribute a partial estimate of the exposure-time surface, while
#' controlling for a model given by 'formula'.
#'
#' @param formula object of class formula, a symbolic description of the fixed
#' effect model to be fitted, e.g. y ~ a + b
#' @param data data frome containing variables used in the formula
#' @param exposure.data numerical matrix of complete exposure data with same
#' length as data
#' @param exposure.splits scalar indicating the number of splits (divided
#' evenly across quantiles of the exposure data) or list with two components:
#' 'type' = 'values' or 'quantiles', and 'split.vals' = a numerical
#' vector indicating the corresponding exposure values or quantiles for splits.
#' Setting exposure.splits equal to 0 will change the model to a distributed lag
#' model, which assumes a linear effect of exposure.
#' @param exposure.se numerical matrix of exposure standard errors with same
#' size as exposure.data or a scalar smoothing factor representing a uniform
#' smoothing factor applied to each exposure measurement
#' @param n.trees integer for number of trees in ensemble
#' @param n.burn integer for length of burn-in
#' @param n.iter integer for number of iterations to run model after burn-in
#' @param n.thin integer thinning factor, i.e. keep every tenth iteration
#' @param family 'gaussian' for continuous response, 'logit' for binomial
#' response with logit link, 'zinb' for zero-inflated negative binomial model
#' @param binomial.size integer type scalar (if all equal, default = 1) or
#' vector defining binomial size for 'logit' family
#' @param tree.params numerical vector of alpha and beta hyperparameters
#' controlling tree depth (see Bayesian CART, 1998), default: alpha = 0.95,
#' beta = 2
#' @param step.prob numerical vector for probability of 1) grow/prune, and
#' 2) change, defaults to (0.25, 0.25) or equal
#' probability of each step for tree updates
#' @param shrinkage int, 1 (default) turn on tree-specific shrinkage priors,
#' 0 turn off
#' @param subset integer vector to analyze only a subset of data and exposures
#' @param verbose true (default) or false: print output
#' @param diagnostics true or false (default) keep model diagnostic such as
#' terminal nodes, acceptance details, etc.
#' @param ... NA
#'
#' @details Model is recommended to be run for at minimum 5000 burn-in iterations
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
                  family = "gaussian",
                  binomial.size = 1,
                  tree.params = c(.95, 2),
                  step.prob = c(.25, .25),
                  shrinkage = 1,
                  subset = NULL,
                  verbose = TRUE,
                  diagnostics = FALSE,
                  spNodes = NULL,
                  areaW = NULL,
                  areaA = NULL,
                  ...)
{
  model <- list()
  options(stringsAsFactors = F)

  # ---- Check inputs ----
  # number of iterations
  if (n.iter < n.thin * 10)
    stop("after thinning you will be left with less than 10 MCMC samples,
         increase the number of iterations!")

  # data for formula
  if (!is.data.frame(data))
    stop("`data` must be a data.frame")

  # exposure data (single exposure)
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
    if (!all(dim(exposure.data) == dim(exposure.se)))
      stop("`exposure.se` and `exposure.data` must have same dimensions")
  }

  # iteration control
  if (all(sapply(list(n.trees, n.burn, n.iter, n.thin),
                 function(i) is.integer(i) & i > 0)))
    stop("n.* must be integer and > 0")

  # response type
  if (!(family %in% c("gaussian", "logit", "zinb")))
    stop("`family` must be one of `gaussian`, `logit`, or `zinb`")

  # binomial size
  model$binomial <- 0
  model$binomialSize <- rep(0, nrow(data))
  if (family == "logit") {
    if (length(binomial.size) == 1)
      binomial.size <- rep(binomial.size, nrow(data))
    if (length(binomial.size) != nrow(data))
      stop("`binomial.size` must be positive integer and same length as data")
    model$binomialSize <- force(binomial.size)
    model$binomial <- 1
  }

  # If not ZINB, set a flag to false (0)
  model$zinb <- 0
  # If ZINB is called, set the flag to true
  if(family == "zinb"){
    model$zinb <- 1
    #model$wTrue = wTrue
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


  # ---- Model control arguments ----
  model$nTrees <- force(n.trees)
  model$nBurn <- force(n.burn)
  model$nIter <- force(n.iter)
  model$nThin <- force(n.thin)
  model$mcmcIter <- force(floor(n.iter / n.thin))
  model$family <- force(family) 
  model$verbose <- force(verbose)
  model$diagnostics <- force(diagnostics)
  model$treePrior <- force(tree.params)
  model$stepProb <- force(c(step.prob[1], step.prob[1], step.prob[2]))
  model$stepProb <- force(model$stepProb / sum(model$stepProb))
  model$shrinkage <- shrinkage

  if (model$verbose)
    cat("Preparing data...\n")



  # ---- Create data subset ----
  if (!is.null(subset)) {
    if (length(subset) > 1 & is.integer(subset) &
        all(subset > 0) & all(subset <= nrow(data))) {
      data <- data[subset,]
      exposure.data <- exposure.data[subset,]
      if (!is.null(exposure.se))
        exposure.se <- exposure.se[subset,]
      if (model$family == "logit")
        model$binomialSize <- model$binomialSize[subset]
    } else {
      stop("invalid subset, must be integers within range of data length")
    }
  }


  # ---- Setup control and response variables ----
  model$formula <- force(as.formula(formula))
  tf <- terms.formula(model$formula, data = data)
  if (!attr(tf, "response"))
    stop("no valid response in formula")
  model$intercept <- force(ifelse(attr(tf, "intercept"), TRUE, FALSE))
  if (length(which(attr(tf, "term.labels") %in% colnames(data))) == 0)
    stop("no valid variables in formula")


  # ---- Exposure splits ----
  model$pExp <- ncol(exposure.data)
  model$timeProb <- force(rep(1 / (model$pExp - 1), model$pExp - 1))
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
      model$X <- force(model$X / model$Xscale)
      model$Tcalc <- force(sapply(1:ncol(model$X),
                                  function(i) rowSums(model$X[, 1:i, drop = F])))

    # splits defined by quantiles of exposure
    } else {
      model$dlFunction <- "tdlnm"
      if (is.list(exposure.splits)) {
        stop("exposure.splits must be a scalar or list with two inputs:
             'type' and 'split.vals'")
      } else {
        model$Xsplits <- force(sort(unique(quantile(model$X,
                                                    (1:(exposure.splits - 1)) /
                                                      exposure.splits))))
        model$nSplits <- force(length(model$Xsplits))
        model$splitProb <- force(rep(1 / model$nSplits, model$nSplits))
      }
    }

  # splits defined by specific values or quantiles
  } else {

    model$dlFunction <- "tdlnm"

    # if exposure.splits entered incorrectly, infer user input and inform
    if (!is.list(exposure.splits)) {
      if (any(exposure.splits > 1 | exposure.splits < 0)) {
        cat("exposure.splits entered as numeric vector, assuming values are
            exposure splitting points\n")
        exposure.splits <- list("type" = "values",
                                "split.vals" = exposure.splits)
      } else {
        cat("exposure.splits entered as numeric vector, assuming values are
            exposure splitting quantiles\n")
        exposure.splits <- list("type" = "quantiles", "split.vals" = exposure.splits)
      }
    }

    # use specific values as splitting points
    if (exposure.splits$type == "values") {
      model$Xsplits <- force(sort(unique(exposure.splits$split.vals)))
      model$Xsplits <- force(model$Xsplits[which(model$Xsplits > min(model$X) &
                                                 model$Xsplits < max(model$X))])

    # use specific quantiles as splitting points
    } else if (exposure.splits$type == "quantiles") {
      if (any(exposure.splits$split.vals > 1 | exposure.splits$split.vals < 0))
        stop("`exposure.splits$split.vals` must be between zero and one if
             using quantiles")
      model$Xsplits <- force(sort(unique(quantile(model$X,
                                                  exposure.splits$split.vals))))
      model$Xsplits <- force(model$Xsplits[which(model$Xsplits > min(model$X) &
                                                 model$Xsplits < max(model$X))])
    } else {
      stop("`exposure.splits$type` must be one of `values` or `quantiles`")
    }

    model$nSplits <- force(length(model$Xsplits))
    if (model$nSplits == 0)
      stop("no exposure splits specified, please check `exposure.splits` input")

    model$splitProb <- force(rep(1 / model$nSplits, model$nSplits))
  }


  # Precalculate counts below each splitting values
  if (length(model$Xsplits) > 0) {
    model$Xscale <- 1
    model$Tcalc <- force(sapply(1:ncol(model$X), function(i) {
      rep(i / model$Xscale, nrow(model$X)) }))
    if (model$smooth) {
      model$Xcalc <- force(sapply(model$Xsplits, function(i) {
        rowSums(pnorm((i - model$X) / model$SE)) }) / model$Xscale)
    } else {
      model$Xcalc <- force(sapply(model$Xsplits, function(i) {
        rowSums(model$X < i) }) / model$Xscale)
    }
  }



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
  rm(QR)



  if (model$family == "gaussian") {
    model$Ymean <- sum(range(model$Y))/2
    model$Yscale <- diff(range(model$Y - model$Ymean))
    model$Y <- force((model$Y - model$Ymean) / model$Yscale)
  } else { # logistic & Gaussian
    model$Yscale <- 1
    model$Ymean <- 0
    model$Y <- force(scale(model$Y, center = 0, scale = 1))
  }
  model$Y <- force(c(model$Y))
  model$Zscale <- attr(model$Z, "scaled:scale")
  model$Zmean <- attr(model$Z, "scaled:center")
  model$Znames <- colnames(model$Z)
  model$Z <- force(matrix(model$Z, nrow(model$Z), ncol(model$Z)))

  # ---- [Spatial structure construction: Adjacency & Diagonal matrix for CAR model] ----
  # default setting for no spatial effect
  model$spatial = FALSE

  #if(!is.null(spNodes) && !is.null(areaA)){
  if(!is.null(spNodes) && !is.null(areaW) && !is.null(areaA)){
    # For now, we just feed the cov-var matrix of phi for the compuation issue with spatial packages in HPC
    model$spNodes1 = spNodes[[1]]
    model$spNodes2 = spNodes[[2]]
    model$areaW = areaW        
    model$areaD = diag(rowSums(areaW))

    model$areaA = force(scaleModelMatrix(areaA))      # Assignment matrix A so that A*phi is n x 1
    model$areaA = force(matrix(model$areaA, nrow(model$areaA), ncol(model$areaA)))
    model$Ascale <- attr(model$areaA, "scaled:scale")
    model$Amean <- attr(model$areaA, "scaled:center")

    model$spN = ncol(areaA)                           # Number of unique areas
    model$spatial = TRUE                              # Set the spatial flag to be true
  }

  # ---- Run model ----
  # model.list <- lapply(ls(envir = model), function(i) model[[i]])
  # names(model.list) <- ls(envir = model)
  out <- tdlnm_Cpp(model)

  # rm(model.list)
  if (verbose)
    cat("\nCompiling results...")

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
  # Gaussian & Logistic fixed effect
  model$gamma <- sapply(1:ncol(model$gamma), function(i) {
  model$gamma[,i] * model$Yscale / model$Zscale[i] })

  # ZINB fixed effect (binary)
  model$b1 <- sapply(1:ncol(model$b1), function(i) {
  model$b1[,i] * model$Yscale / model$Zscale[i] })

  # ZINB fixed effect (NegBin)
  model$b2 <- sapply(1:ncol(model$b2), function(i) {
  model$b2[,i] * model$Yscale / model$Zscale[i] })

  if (model$intercept) {
    model$gamma[,1] <- model$gamma[,1] + model$Ymean
    if (ncol(model$Z) > 1)
      model$gamma[,1] <- model$gamma[,1] - model$gamma[,-1, drop=FALSE] %*% model$Zmean[-1]

    model$b1[,1] <- model$b1[,1] + model$Ymean
    if (ncol(model$Z) > 1)
      model$b1[,1] <- model$b1[,1] - model$b1[,-1] %*% model$Zmean[-1]

    model$b2[,1] <- model$b2[,1] + model$Ymean
    if (ncol(model$Z) > 1)
      model$b2[,1] <- model$b2[,1] - model$b2[,-1] %*% model$Zmean[-1]
  }

  colnames(model$gamma) <- model$Znames
  colnames(model$b1) <- model$Znames
  colnames(model$b2) <- model$Znames


  # rescale DLM estimates
  model$DLM <- as.data.frame(model$DLM)
  colnames(model$DLM) <- c("Iter", "Tree", "xmin", "xmax", "tmin", "tmax",
                           "est", "kappa")
  model$DLM$est <- model$DLM$est * model$Yscale / model$Xscale
  model$DLM$xmin <- sapply(model$DLM$xmin, function(i) {
    if (i == 0) -Inf
    else model$Xsplits[i]
  })
  model$DLM$xmax <- sapply(model$DLM$xmax, function(i) {
    if (i == (length(model$Xsplits) + 1)) Inf
    else model$Xsplits[i]
  })

  # Remove model and exposure data
  model$X <- NULL
  model$Tcalc <- NULL
  model$Xcalc <- NULL
  model$Z <- NULL

  # Change env to list
  model.out <- lapply(names(model), function(i) model[[i]])
  names(model.out) <- names(model)
  rm(model)
  gc()

  class(model.out) <- ifelse(model.out$dlFunction == "tdlnm", "tdlnm", "tdlm")
  return(model.out)
}
