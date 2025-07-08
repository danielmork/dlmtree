#' dlmtree
#'
#' @title Fit tree structured distributed lag models
#' @description 
#' The 'dlmtree' function accommodates various response variable types, including continuous, binary, and zero-inflated count values. 
#' The function is designed to handle both single exposure and exposure mixtures. For a single exposure, users are offered options to model non-linear effects (tdlnm), 
#' linear effects (tdlm), or heterogeneous subgroup/individualized effects (hdlm). In the case of exposure mixtures, the function supports 
#' lagged interactions (tdlmm), and heterogeneous subgroup/individualized effects (hdlmm) allowing for a comprehensive exploration of mixture exposure heterogeneity. 
#' Additionally, users can fine-tune parameters to impose effect shrinkage and perform exposure selection, enhancing the adaptability and precision of the modeling process.
#' For more detailed documentation, visit: [dlmtree website](https://danielmork.github.io/dlmtree/).
#'
#' @param formula object of class formula, a symbolic description of the fixed
#' effect model to be fitted, e.g. y ~ a + b.
#' @param data data frame containing variables used in the formula.
#' @param exposure.data numerical matrix of exposure data with same length as data, for a mixture setting (tdlmm, hdlmm): 
#' named list containing equally sized numerical matrices of exposure data having same length as data.
#' @param dlm.type dlm model specification: "linear" (default), "nonlinear", "monotone".
#' @param family 'gaussian' for continuous response, 'logit' for binomial, 'zinb' for zero-inflated negative binomial.
#' @param mixture flag for mixture, set to TRUE for tdlmm and hdlmm (default: FALSE).
#' @param het flag for heterogeneity, set to TRUE for hdlm and hdlmm (default: FALSE).
#' @param control.mcmc list of MCMC control parameters. This is passed to \link{dlmtree.control.mcmc}.
#' @param control.hyper list of hyperparameter control parameters. This is passed to \link{dlmtree.control.hyper}
#' @param control.family list of family control parameters. This is passed to \link{dlmtree.control.family}
#' @param control.tdlnm list of TDLNM control parameters. This is passed to \link{dlmtree.control.tdlnm}
#' @param control.het list of control parameters for heterogeneous models. This is passed to \link{dlmtree.control.het}
#' @param control.mix list of mixture control parameters. This is passed to \link{dlmtree.control.mix}
#' @param control.monotone list of control parameters for monotone model. This is passed to \link{dlmtree.control.monotone}
#' @param control.diagnose list of control parameters for diagnostics. This is passed to \link{dlmtree.control.diagnose}
#'
#' @details Model is recommended to be run for at minimum 5000 burn-in
#' iterations followed by 15000 sampling iterations with a thinning factor of 5.
#' Convergence can be checked by re-running the model and validating consistency
#' of results. Examples are provided below for the syntax for running different types of models.
#' For more examples, visit: [dlmtree website](https://danielmork.github.io/dlmtree/).
#' 
#' @md
#' @example inst/examples/dlmtree_example.R
#'
#' @returns object of one of the classes: tdlm, tdlmm, tdlnm, hdlm, hdlmm, monotone
#' @export
#'
dlmtree <- function(formula,
                    data,
                    exposure.data,
                    dlm.type = "linear",
                    family = "gaussian",
                    mixture = FALSE,
                    het = FALSE,                
                    # Control parameters
                    control.mcmc = list(),
                    control.hyper = list(),
                    control.family = list(),
                    control.tdlnm = list(),
                    control.het = list(),
                    control.mix = list(),
                    control.monotone = list(),
                    control.diagnose = list())
{
  # *** Extract control parameters ***
  control.mcmc     <- do.call("dlmtree.control.mcmc", control.mcmc)
  control.hyper    <- do.call("dlmtree.control.hyper", control.hyper)
  control.family   <- do.call("dlmtree.control.family", control.family)
  control.tdlnm    <- do.call("dlmtree.control.tdlnm", control.tdlnm)
  control.het      <- do.call("dlmtree.control.het", control.het)
  control.mix      <- do.call("dlmtree.control.mix", control.mix)
  control.monotone <- do.call("dlmtree.control.monotone", control.monotone)
  control.diagnose <- do.call("dlmtree.control.diagnose", control.diagnose)
  
  
  # MCMC control parameters
  n.trees <- control.mcmc$n.trees
  n.burn  <- control.mcmc$n.burn
  n.iter  <- control.mcmc$n.iter
  n.thin  <- control.mcmc$n.thin
  
  # Hyperparameter control parameters
  shrinkage         <- control.hyper$shrinkage
  dlmtree.params    <- control.hyper$params
  dlmtree.step.prob <- control.hyper$step.prob
  
  # Family control parameters
  binomial.size <- control.family$binomial.size
  formula.zi    <- control.family$formula.zi
  
  # TDLNM control parameters
  tdlnm.exposure.splits <- control.tdlnm$exposure.splits
  tdlnm.time.split.prob <- control.tdlnm$time.split.prob
  tdlnm.exposure.se     <- control.tdlnm$exposure.se
  
  # Heterogeneous model control parameters
  hdlm.modifiers          <- control.het$modifiers
  hdlm.modifier.splits    <- control.het$modifier.splits
  hdlm.modtree.params     <- control.het$modtree.params
  hdlm.modtree.step.prob  <- control.het$modtree.step.prob
  hdlm.dlmtree.type       <- control.het$dlmtree.type
  hdlm.selection.prior    <- control.het$selection.prior
  
  # Mixture model control parameters
  mixture.interactions  <- control.mix$interactions
  mixture.prior         <- control.mix$sparsity.prior
  
  # Monotone model control parameters
  monotone.gamma0           <- control.monotone$gamma0
  monotone.sigma            <- control.monotone$sigma
  monotone.tree.time.params <- control.monotone$tree.time.params
  monotone.tree.exp.params  <- control.monotone$tree.exp.params
  monotone.time.kappa       <- control.monotone$time.kappa
  
  # Diagnostic control parameters
  subset          <- control.diagnose$subset
  lowmem          <- control.diagnose$lowmem
  verbose         <- control.diagnose$verbose
  save.data       <- control.diagnose$save.data
  diagnostics     <- control.diagnose$diagnostics
  initial.params  <- control.diagnose$initial.params


  # *** Initialize model ***
  model <- list()
  options(stringsAsFactors = FALSE)
  piecewise.linear = FALSE

  
  # *** Check MCMC inputs ***
  if (all(sapply(list(n.trees, n.burn, n.iter, n.thin), function(i) is.integer(i) & i > 0))) {
    stop("n.* must be integer and > 0")
  }

  if (n.iter < n.thin * 10) {
    stop("After thinning, you will be left with less than 10 MCMC samples,",
          " increase the number of iterations!")
  }
  
  
  # *** Check model specification ***
  # family
  if (!(family %in% c("gaussian", "logit", "zinb"))) {
    stop("`family` must be one of `gaussian`, `logit`, or 'zinb'")
  }

  # dlm.type
  if (!(dlm.type %in% c("linear", "nonlinear", "monotone"))) {
    stop("`dlm.type` must be one of `linear`, `nonlinear`, `monotone`")
  }

  # check shared / nested
  if (!(hdlm.dlmtree.type %in% c("shared", "nested"))){
    stop("`hdlm.dlmtree.type` must be one of `shared`, `nested`")
  }

  # Stop for unavailable models
  if (het) { # HDLM & HDLMM
    if (family %in% c("logit", "zinb")) {
      stop("'logit' or 'zinb' are unavailable for heterogeneous models. Set family to 'gaussian'.")
    }

    if (dlm.type %in% c("nonlinear", "monotone")) {
      stop("Non-linear or monotone DLM types are unavailable for heterogeneous models. Set dlm.type to 'linear'.")
    }
  } else { # TDLM, TDLNM, TDLMM
    if (dlm.type != "linear"){
      if (family == "zinb") {
        stop("'zinb' is unavailable for non-linear or monotone models. Set family to 'gaussian'.")
      }
    }

    if(dlm.type %in% c("nonlinear", "monotone") & sum(tdlnm.exposure.splits) == 0){
      stop("Exposure split cannot be set as zero for non-linear or monotone models. 
            Set dlm.type to 'linear', which will result in running the TDLM model")
    }
  }


  # *** Check data inputs ***
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame")
  }

  # Single exposure models (TDLM, TDLNM, Monotone)
  if (!mixture) {
    
    if (!is.numeric(exposure.data)) {
      stop("`exposure.data` must be a numeric matrix for single exposure models")
    }
      
    if (nrow(data) != nrow(exposure.data)) {
      stop("`data` and `exposure.data` must have same number of observations")
    }
      
    model$pExp <- ncol(exposure.data) # lag

    # TDLM 
    if (dlm.type == "linear") {
      tdlnm.exposure.splits <- 0
    }

    model$monotone <- ifelse(dlm.type == "monotone", TRUE, FALSE)

    # Monotone default set up
    if(is.null(monotone.gamma0)){
      monotone.gamma0 <- rep(0, ncol(exposure.data))
    }

    if(is.null(monotone.sigma)){
      monotone.sigma <- diag(ncol(exposure.data)) * 1.502^2
    }

    # Exposure smoothing/error (TDLNM, monotone)
    # Set default
    if(is.null(tdlnm.exposure.se)){
      tdlnm.exposure.se <- sd(exposure.data)/2
    }

    if(is.null(tdlnm.time.split.prob)){
      tdlnm.time.split.prob <- rep(1 / (model$pExp - 1), model$pExp - 1)
    }

    # check the dimension
    if (!is.numeric(tdlnm.exposure.se)) {
      stop("`tdlnm.exposure.se` must be a scalar or numeric matrix")
    }

    if (length(tdlnm.exposure.se) == 1) {
      tdlnm.exposure.se <- matrix(tdlnm.exposure.se, nrow(exposure.data), ncol(exposure.data))
    }
    
    if (!all(dim(exposure.data) == dim(tdlnm.exposure.se))) {
      stop("`tdlnm.exposure.se` and `exposure.data` must have same dimensions")
    }

  } else { # HDLMM, TDLMM
    exp.centered <- FALSE # Check exposure centering for TDLMM

    # TDLNM (no mixture) 
    if (dlm.type == "nonlinear") {
      stop("The mixture setting is not applicable to `non-linear` model. Set mixture to FALSE.")
    }

    if (!is.list(exposure.data)) {
      stop("`exposure.data` must be a list of named exposures for mixture models")
    }

    model$nExp <- length(exposure.data)       # Number of exposures
    model$expNames <- names(exposure.data)    # Exposure names

    if (is.null(model$expNames) || length(unique(model$expNames)) != model$nExp || any(model$expNames == "")) {
      stop("`exposure.data` must be a named list with unique, non-empty names")
    }

    model$pExp <- ncol(exposure.data[[1]])    # total lag

    for(i in 1:length(exposure.data)) {
      if (!is.numeric(exposure.data[[i]])) {
        stop("each exposure in list `exposure.data` must be a numeric matrix")
      }

      if (nrow(data) != nrow(exposure.data[[i]])) {
        stop("`data` and `exposure.data` must have same number of rows")
      }
        
      if (ncol(exposure.data[[i]]) != model$pExp) {
        stop("each exposure in `exposure.data` must have the same number of time observations")
      }

      if (any(is.na(exposure.data[[i]]))) {
        stop("missing values in exposure data")
      }

      # Check exposure centering for TDLMM
      if (het == FALSE & mean(exposure.data[[i]]) < 1e-10){
        exp.centered <- TRUE
      }
    }

    if(exp.centered){
      message("Caution: At least one of the exposure data may have been centered. This could result in an inaccurate estimate of marginal exposure effect when using TDLMM \n")
    }
  }

  # *** Check response type ***
  # Binary flag
  model$binomial <- 0
  model$binomialSize <- rep(0, nrow(data))
  if (family == "logit") {
    if (length(binomial.size) == 1)
      binomial.size <- rep(binomial.size, nrow(data))
    if (length(binomial.size) != nrow(data))
      stop("`binomial.size` must be positive integer and same length as data")
    model$binomialSize <- binomial.size
    model$binomial <- 1
  }

  # Zero-inflated negative binomial flag
  model$zinb <- 0
  if (family == "zinb") {
    model$zinb    <- 1
    model$sigma2  <- 1
  }
  
  # Mixture interactions
  if (!(mixture.interactions %in% c("none", "noself", "all"))) {
    stop("`mixture.interactions must be one of `none`, `noself`, `all`")
  }

  if (mixture.interactions %in% c("marginal", "none")) {
    model$interaction <- 0     
    model$nMix        <- 0   
  } else if (mixture.interactions == "noself") { 
    model$interaction <- 1  
    model$nMix        <- model$nExp * (model$nExp - 1) / 2 
  } else { 
    model$interaction <- 2   
    model$nMix        <- model$nExp * (model$nExp + 1) / 2 
  }

  # *** Check tree parameters ***
  # Prior splitting probability
  if (length(hdlm.modtree.params) != 2 || length(dlmtree.params) != 2) {
    stop("tree.params.* must have length 2")
  } else {
    if (hdlm.modtree.params[1] > 1 || hdlm.modtree.params[1] < 0 || dlmtree.params[1] > 1 || dlmtree.params[1] < 0) {
          stop("tree.params.*[1] must be between 0 and 1")
        }

    if (hdlm.modtree.params[2] < 0 || dlmtree.params[2] < 0) {
      stop("tree.params.*[2] must be greater than 0")
    }
  }

  # Step probabilities
  if (any(hdlm.modtree.step.prob < 0) || any(hdlm.modtree.step.prob > 1) || any(dlmtree.step.prob < 0) || any(dlmtree.step.prob > 1)) {
    stop("step.prob.* components must be between zero and one")
  }

  
  # *** Model control arguments ***
  # MCMC
  model$nTrees    <- n.trees
  model$nBurn     <- n.burn
  model$nIter     <- n.iter
  model$nThin     <- n.thin
  model$mcmcIter  <- floor(n.iter / n.thin)
  
  # Model specification
  model$family    <- family
  model$mixture   <- mixture
  model$het       <- het

  # dlmtree
  model$treePriorTDLM <- dlmtree.params
  model$stepProbTDLM  <- prop.table(c(dlmtree.step.prob[1], dlmtree.step.prob[1], dlmtree.step.prob[2]))
  if(mixture){
    model$stepProbTDLM  <- prop.table(c(dlmtree.step.prob[1], dlmtree.step.prob[1], dlmtree.step.prob[2], dlmtree.step.prob[2]))
  }
  model$timeSplits0 <- tdlnm.time.split.prob

  # modtree
  model$treePriorMod  <- hdlm.modtree.params
  model$stepProbMod   <- prop.table(c(hdlm.modtree.step.prob[1], hdlm.modtree.step.prob[2], hdlm.modtree.step.prob[3], 1 - sum(hdlm.modtree.step.prob)))

  # Monotone model priors            
  if (!mixture) { # if statement to avoid a warning when running TDLMM
    model$shape     <- ifelse(mean(tdlnm.exposure.splits) == 0, "Linear",
                              ifelse(mean(tdlnm.exposure.se) != 0, "Smooth",
                                     "Step Function"))
  }
  
  if (dlm.type == "monotone") {
    model$zirtGamma0      <-  monotone.gamma0
    model$zirtSigma       <-  monotone.sigma
    model$treePriorTime   <-  monotone.tree.time.params
    model$treePriorExp    <-  monotone.tree.exp.params
    model$timeKappa       <-  ifelse(is.null(monotone.time.kappa), 1.0, monotone.time.kappa)
    model$updateTimeKappa <-  ifelse(is.null(monotone.time.kappa), TRUE, FALSE)
  }

  # Modifier prior for HDLM, HDLMM
  if (hdlm.selection.prior < 0.5 | hdlm.selection.prior > 1){
    stop("'hdlm.selection.prior' must be between 0.5 and 1, inclusive")
  }
  model$zeta <- hdlm.selection.prior

  # Mixture parameters
  if (!(shrinkage %in% c("all", "trees", "exposures", "none"))) {
    stop("`shrinkage` must be one of `all`, `trees`, `exposures`, or `none`")
  }

  model$mixPrior  <- mixture.prior # Exposure selection sparsity
  model$shrinkage <- ifelse(dlm.type == "monotone", FALSE,
                      switch(shrinkage, "all" = 3, "trees" = 2, "exposures" = 1, "none" = 0))

  # Diagnostics
  model$lowmem      <- lowmem
  model$verbose     <- verbose
  model$diagnostics <- diagnostics
  model$debug       <- FALSE
  
  if (verbose) {
    cat("Preparing data...\n")
  }

  # *** Create data subset ***
  if (!is.null(subset)) {
    if (!mixture) {
      if (!is.null(subset)) {
        if (length(subset) > 1 & is.integer(subset) & all(subset > 0) & all(subset <= nrow(data))) {
          data <- data[subset,]
          exposure.data <- exposure.data[subset,]

          if (model$shape != "Linear")
            tdlnm.exposure.se <- tdlnm.exposure.se[subset,]

          if (model$family == "logit")
            model$binomialSize <- model$binomialSize[subset]

        } else {
          stop("invalid subset, must be integers within range of data length")
        }
      }

    } else {
      if (!is.integer(subset) || any(subset < 1) || any(subset > nrow(data))) {
        stop("invalid subset, must be integers within range of data length")
      }
        
      data <- data[subset,]
      exposure.data <- lapply(exposure.data, function(i) i[subset,])
    }
  }

  # *** Setup fixed effect model ***
  if (is.null(formula.zi)) {
    formula.zi <- formula
  }

  model$formula <- as.formula(formula)
  model$formula.zi  <- as.formula(formula.zi)

  tf   <- terms.formula(model$formula, data = data)
  tf.zi    <- terms.formula(model$formula.zi, data = data) 

  # Sanity check for response for variables and the names
  if (!attr(tf, "response") & !attr(tf.zi, "response")){
    stop("no valid response in formula")
  }

  model$intercept <- ifelse(attr(tf, "intercept"), TRUE, FALSE)
  model$intercept.zi <- ifelse(attr(tf.zi, "intercept"), TRUE, FALSE)

  # Sanity check for variables and the names
  if (length(which(attr(tf, "term.labels") %in% colnames(data))) == 0 & 
      length(which(attr(tf.zi, "term.labels") %in% colnames(data))) == 0 & 
      !model$intercept & !model$intercept.zi) 
    stop("no valid variables in formula")

  # *** Exposure splits ***
  model$splitProb <- as.double(c())
  model$timeProb  <- rep(1 / (model$pExp - 1), model$pExp - 1)
  model$nSplits   <- 0

  # Model processing
  if (het) {
    if (mixture) { # HDLMM
      model$class <- "hdlmm"
      model$X     <- list() 
      for(i in 1:model$nExp) { # For each exposure,
        model$X[[i]]        <- list(Xscale = sd(exposure.data[[i]]), X = exposure.data[[i]])
        model$X[[i]]$X      <- model$X[[i]]$X / model$X[[i]]$Xscale
        model$X[[i]]$Xrange <- range(model$X[[i]]$X)
        model$X[[i]]$Xquant <- quantile(model$X[[i]]$X, 0:100/100) *  model$X[[i]]$Xscale
        model$X[[i]]$intX   <- mean(model$X[[i]]$X)
        model$X[[i]]$Tcalc  <- sapply(1:ncol(model$X[[i]]$X), function(j) {rowSums(model$X[[i]]$X[, 1:j, drop = FALSE])}) 
      }

      names(model$X)        <- model$expNames
      model$expProb         <- rep(1/length(model$X), length(model$X))
    } else { # HDLM
      model$class   <- "hdlm"
      model$X       <- exposure.data
      model$Xrange  <- range(model$X)
      model$Xscale  <- sd(model$X)
      model$X       <- model$X / model$Xscale
      model$Tcalc   <- sapply(1:ncol(model$X), function(i) rowSums(model$X[, 1:i, drop = FALSE]))
    }
  } else {
    if (mixture) { # TDLMM
      model$class <- "tdlmm"
      model$X     <- list() 
      for(i in 1:model$nExp) { # For each exposure,
        model$X[[i]]        <- list(Xscale = sd(exposure.data[[i]]), X = exposure.data[[i]])
        model$X[[i]]$X      <- model$X[[i]]$X / model$X[[i]]$Xscale
        model$X[[i]]$Xrange <- range(model$X[[i]]$X)
        model$X[[i]]$Xquant <- quantile(model$X[[i]]$X, 0:100/100) *  model$X[[i]]$Xscale
        model$X[[i]]$intX   <- mean(model$X[[i]]$X)
        model$X[[i]]$Tcalc  <- sapply(1:ncol(model$X[[i]]$X), function(j) {rowSums(model$X[[i]]$X[, 1:j, drop = FALSE])}) 
      }
      names(model$X)        <- model$expNames
      model$expProb         <- rep(1/length(model$X), length(model$X))
    } else {
      model$X       <- exposure.data
      model$Xrange  <- range(exposure.data)

      if (mean(tdlnm.exposure.se) == 0) {
        model$smooth <- FALSE
        model$SE <- matrix(0.0, 0, 0)
      } else {
        model$smooth <- TRUE
        model$SE <- tdlnm.exposure.se
      }

      if (length(tdlnm.exposure.splits) == 1) {
        # TDLM
        if (tdlnm.exposure.splits == 0) {
          model$class     <- "tdlm"
          model$splitProb <- as.double(c())
          model$Xsplits   <- as.double(c())
          model$nSplits   <- 0
          model$Xscale    <- sd(model$X)
          model$X         <- model$X / model$Xscale
          model$Tcalc     <- sapply(1:ncol(model$X),
                                      function(i) rowSums(model$X[, 1:i, drop = FALSE]))

        # TDLNM: Splits defined by quantiles of exposure
        } else {
          model$class <- ifelse(model$monotone, "monotone", "tdlnm")
          
          if (is.list(tdlnm.exposure.splits)) {
            stop("tdlnm.exposure.splits must be a scalar or list with two inputs: 'type' and 'split.vals'")
          } else {
            model$Xsplits   <- sort(unique(quantile(model$X,
                                                        (1:(tdlnm.exposure.splits - 1)) /
                                                          tdlnm.exposure.splits)))
            model$nSplits   <- length(model$Xsplits)
            model$splitProb <- rep(1 / model$nSplits, model$nSplits)
          }
        }

      # TDLNM: Splits defined by specific values or quantiles
      } else {
        model$class <- ifelse(model$monotone, "monotone", "tdlnm")

        # if tdlnm.exposure.splits entered incorrectly, infer user input and inform
        if (!is.list(tdlnm.exposure.splits)) {
          if (any(tdlnm.exposure.splits > 1 | tdlnm.exposure.splits < 0)) {
            if (verbose)
              cat("tdlnm.exposure.splits entered as numeric vector, assuming values are exposure splitting points\n")
            tdlnm.exposure.splits <- list("type" = "values",
                                          "split.vals" = tdlnm.exposure.splits)
          } else {
            if (verbose)
              cat("tdlnm.exposure.splits entered as numeric vector, assuming values are exposure splitting quantiles\n")
            tdlnm.exposure.splits <- list("type" = "quantiles", "split.vals" = tdlnm.exposure.splits)
          }
        }

        # use specific values as splitting points
        if (tdlnm.exposure.splits$type == "values") {
          model$Xsplits <- sort(unique(tdlnm.exposure.splits$split.vals))
          model$Xsplits <- model$Xsplits[which(model$Xsplits > min(model$X) &
                                                    model$Xsplits < max(model$X))]

        # use specific quantiles as splitting points
        } else if (tdlnm.exposure.splits$type == "quantiles") {
          if (any(tdlnm.exposure.splits$split.vals > 1 | tdlnm.exposure.splits$split.vals < 0))
            stop("`tdlnm.exposure.splits$split.vals` must be between zero and one if using quantiles")
          model$Xsplits <- sort(unique(quantile(model$X,
                                                      tdlnm.exposure.splits$split.vals)))
          model$Xsplits <- model$Xsplits[which(model$Xsplits > min(model$X) &
                                                    model$Xsplits < max(model$X))]
        } else {
          stop("`tdlnm.exposure.splits$type` must be one of `values` or `quantiles`")
        }

        model$nSplits <- length(model$Xsplits)
        if (model$nSplits == 0)
          stop("no exposure splits specified, please check `tdlnm.exposure.splits` input")

        model$splitProb <- rep(1 / model$nSplits, model$nSplits)
        
        # memory warning
        if (prod(dim(model$X)) * model$nSplits * 8 > 1024^3 & model$verbose)
          warning(paste0("Model run will require at least ", 
                        round(prod(dim(model$X)) * model$nSplits * 8 / 1024^3, 1),
                        " GB of memory. Use `lowmem = TRUE` option to reduce memory usage."))
      }
    }
  }

  # Precalculate counts below each splitting values
  if (length(model$Xsplits) > 0) {
    model$Xscale <- 1
    model$Tcalc <- sapply(1:ncol(model$X), function(i) {
      rep(i / model$Xscale, nrow(model$X)) })
    if (model$smooth) {
      model$Xcalc <- sapply(model$Xsplits, function(i) {
        rowSums(pnorm((i - model$X) / model$SE)) }) / model$Xscale
    } else {
      model$Xcalc <- sapply(model$Xsplits, function(i) {
        rowSums(model$X < i) }) / model$Xscale
    }
  }
  
  # *** Setup modifier variables (for HDLM, HDLMM) ***
  if (het) {
    model$modNames <- hdlm.modifiers 

    if (length(model$modNames) == 1) {
      if (model$modNames == "all") {
        # Take out the unused columns and extract column names
        model$modNames <- attr(tf, "term.labels")
      }
    }

    model$Mo <- lapply(model$modNames, function(m) {
        if (!(m %in% colnames(data))) {
            stop("one or more `hdlm.modifiers` specified is not a column of `data`")
        }
        data[[m]] # Extract data with a string exposure name
      }
    )

    model$MoUnique <- lapply(model$modNames, function(m) {
                                if (!is.numeric(data[[m]])){
                                  sort(unique(data[[m]]))
                                } else {
                                  range(data[[m]])
                                }
                              }
                            )
    names(model$Mo) <- names(model$MoUnique) <- model$modNames
    model$pM        <- length(model$Mo)

    # Use quantiles of modifiers as possible splits
    model$modIsNum        <- sapply(model$Mo, function(i) (is.numeric(i) || is.logical(i)))
    model$modSplitValRef  <- list()
    model$modSplitValIdx  <- list()
    model$modSplitIdx     <- list()
    model$fullIdx         <- 0:(nrow(data) - 1)

    for(i in 1:model$pM) {
      if (model$modIsNum[i]) {
        if (length(unique(model$Mo[[i]])) < hdlm.modifier.splits) {
          uniqueVals                <- sort(unique(model$Mo[[i]]))
          model$modSplitValRef[[i]] <- rowMeans(cbind(uniqueVals[-length(uniqueVals)], uniqueVals[-1]))
        } else {
          uniqueVals                <- sort(unique(quantile(model$Mo[[i]], 1:hdlm.modifier.splits/(hdlm.modifier.splits + 1))))
          model$modSplitValRef[[i]] <- uniqueVals[which(uniqueVals > min(model$Mo[[i]]) & uniqueVals < max(model$Mo[[i]]))]
        }
        model$modSplitValIdx[[i]]   <- 0:(length(model$modSplitValRef[[i]]) - 1)
        model$modSplitIdx[[i]]      <- lapply(model$modSplitValRef[[i]], function(j) which(model$Mo[[i]] < j) - 1)
      } else {
        if (length(unique(model$Mo[[i]])) != length(unique(model$Mo[[i]]))) {
          warning("one or more modifier categories not present in data subset")
        }
        model$modSplitValRef[[i]]   <- sort(unique(model$Mo[[i]]))
        model$modSplitValIdx[[i]]   <- 0:(length(model$modSplitValRef[[i]]) - 1)
        model$modSplitIdx[[i]]      <- lapply(model$modSplitValRef[[i]], function(j) which(model$Mo[[i]] == j) - 1)
      }
    }
  }

  # *** Scale data and setup exposures ***
  data <- droplevels(data)
  mf <- model.frame(tf, data = data, drop.unused.levels = TRUE, na.action = NULL)
  mf.zi <- model.frame(tf.zi, data = data, drop.unused.levels = TRUE, na.action = NULL)

  if (any(is.na(mf)) & any(is.na(mf.zi)))
    stop("missing values in model data, use `complete.cases()` to subset data")

  # Response
  model$Y <- model.response(mf)

  # Check response & model specification
  # Binary response
  if (all(model$Y %in% c(0, 1))) {
    if (family %in% c("gaussian", "zinb")) { # Correct response with wrong model
      warning("The response variable contains only 0s and 1s. Consider using a binomial model (family = `logit`).")
    }
  } else { # Wrong response but correct model
    if (family == "logit") {    
      stop("The response variable contains values other than 0s and 1s, which is not applicable to `logit` model.")
    }
  }

  # Count response
  if (all(sapply(model$Y, function(x) (x %% 1 == 0) & (x >= 0)))) {
    if (family == "gaussian") { # Correct response with wrong model
      if(!all(model$Y %in% c(0, 1))){
        warning("The response variable contains only non-negative integers, 
                Consider using a zero-inflated negative binomial model (family = `zinb`).")
      }
    }
    if (family == "logit" & !all(model$Y %in% c(0, 1))) {
      stop("The response variable contains non-negative integers other than 0s and 1s, 
              Consider using a zero-inflated negative binomial model (family = `zinb`).")
    }
  } else { # Wrong response but correct model
    if (family == "zinb") {
      stop("The response variable contains negative values or non-integer positive values, which is not applicable to `zinb` model.")
    }
  }

  # Covariates for main model
  model$Z <- model.matrix(model$formula, data = mf)
  QR <- qr(crossprod(model$Z))
  model$Z <- model$Z[,sort(QR$pivot[seq_len(QR$rank)])]
  model$droppedCovar <- colnames(model$Z)[QR$pivot[-seq_len(QR$rank)]]
  model$Z <- scaleModelMatrix(model$Z)
  if (length(model$droppedCovar) > 0 & model$verbose) {
    warning("variables {", paste0(model$droppedCovar, collapse = ", "), "} dropped due to perfect collinearity\n")
  }
  rm(QR)

  # Covariates for ZI model
  model$Z.zi <- model.matrix(model$formula.zi, data = mf.zi)
  QR.zi <- qr(crossprod(model$Z.zi)) 
  model$Z.zi <- model$Z.zi[,sort(QR.zi$pivot[seq_len(QR.zi$rank)])]
  model$droppedCovar.zi <- colnames(model$Z.zi)[QR.zi$pivot[-seq_len(QR.zi$rank)]]
  model$Z.zi <- scaleModelMatrix(model$Z.zi)
  if (length(model$droppedCovar.zi) > 0 & model$verbose & family == "zinb") {
    warning("variables {", paste0(model$droppedCovar.zi, collapse = ", "), "} dropped due to perfect collinearity\n")
  }
  rm(QR.zi)

  # *** Scale data ***
  if (model$family == "gaussian") {
    model$Ymean   <- sum(range(model$Y))/2
    model$Yscale  <- diff(range(model$Y - model$Ymean))
    #model$Yscale  <- sd(model$Y - model$Ymean)
    model$Y       <- (model$Y - model$Ymean) / model$Yscale
  } else {
    model$Yscale  <- 1
    model$Ymean   <- 0
    model$Y       <- scale(model$Y, center = 0, scale = 1)
  }

  # Store the processed values
  model$Y <- c(model$Y)

  model$Zscale <- attr(model$Z, "scaled:scale")
  model$Zmean <- attr(model$Z, "scaled:center")
  model$Znames <- colnames(model$Z)
  model$Z <- matrix(model$Z, nrow(model$Z), ncol(model$Z))

  model$Zscale.zi <- attr(model$Z.zi, "scaled:scale")
  model$Zmean.zi <- attr(model$Z.zi, "scaled:center")
  model$Znames.zi <- colnames(model$Z.zi)
  model$Z.zi <- matrix(model$Z.zi, nrow(model$Z.zi), ncol(model$Z.zi))

  # TDLMM
  model$initParams <- rep(0, ncol(model$Z))
  if (!is.null(initial.params)) {
    names(model$initParams) <- colnames(model$Z)
    if (sum(names(initial.params) %in% colnames(model$Z)) > 0) {
      na <- names(initial.params[ # get matching names 
        which(names(initial.params) %in% colnames(model$Z))])
      model$initParams[na] <- initial.params[na]
    }
  }
  
  
  # *** Run model ***
  if (verbose) {
    if(model$class == "monotone"){
      cat(paste0("\nRunning monotone-TDLNM:\n"))
    } else if (model$class == "hdlm"){
      cat(paste0("\nRunning ", hdlm.dlmtree.type, " ", toupper(model$class), ":\n"))
    } else {
      cat(paste0("\nRunning ", toupper(model$class), ":\n"))
    }
  }

  out <- switch(model$class,
                "tdlm"  = tdlnm_Cpp(model),
                "tdlmm" = tdlmm_Cpp(model),
                "hdlm"  = switch(hdlm.dlmtree.type, 
                                 "shared" = dlmtreeHDLMGaussian(model), 
                                 "nested" = dlmtreeTDLM_cpp(model)),
                "hdlmm" = dlmtreeHDLMMGaussian(model),
                "tdlnm" = tdlnm_Cpp(model),
                "monotone" = monotdlnm_Cpp(model))

  if (verbose) {
    cat("\nCompiling results...\n")
  }

  for(n in names(out)) {
    model[[n]] <- out[[n]]
  }

  
  # *** Prepare output ***
  model$Y       <- model$Y * model$Yscale + model$Ymean  
  model$fhat    <- model$fhat * model$Yscale    
  model$sigma2  <- model$sigma2 * (model$Yscale^2)     

  # Shrinkage
  if(shrinkage > 1){
    colnames(model$tau) <- paste0("tau", 1:n.trees)
  }

  # Coefficients
  # Unscale fixed effect estimates
  if (family == "zinb") {
    # ZINB fixed effect (binary)
    model$b1 <- sapply(1:ncol(model$b1), function(i) {
      model$b1[,i] * model$Yscale / model$Zscale.zi[i] })

    # ZINB fixed effect (NB)
    model$b2 <- sapply(1:ncol(model$b2), function(i) {
      model$b2[,i] * model$Yscale / model$Zscale[i] })

    # If intercept:
    if (model$intercept) {
      model$b2[,1] <- model$b2[,1] + model$Ymean

      if (ncol(model$Z) > 1) {
        model$b2[,1] <- model$b2[,1] - model$b2[,-1,drop=FALSE] %*% model$Zmean[-1]
      }
    }

    if (model$intercept.zi) {
      model$b1[,1] <- model$b1[,1] + model$Ymean

      if (ncol(model$Z.zi) > 1) {
        model$b1[,1] <- model$b1[,1] - model$b1[,-1,drop=FALSE] %*% model$Zmean.zi[-1]
      }
    }

    colnames(model$b1) <- model$Znames.zi
    colnames(model$b2) <- model$Znames

  } else { # Gaussian & Logistic fixed effect
    model$gamma <- sapply(1:ncol(model$gamma), function(i) {
      model$gamma[,i] * model$Yscale / model$Zscale[i] })

    if (model$intercept) {
      model$gamma[,1] <- model$gamma[,1] + model$Ymean

      if (ncol(model$Z) > 1) {
        model$gamma[,1] <- model$gamma[,1] - model$gamma[,-1,drop=FALSE] %*% model$Zmean[-1]
      }
    }

    colnames(model$gamma) <- model$Znames
  }

  # Tree information data frames for models
  if (model$class == "tdlmm") {
    # rescale DLM and Mixture estimates
    model$TreeStructs           <- as.data.frame(model$TreeStructs)
    colnames(model$TreeStructs) <- c("Iter", "Tree", "TreePair", "exp", "tmin", "tmax", "est", "kappa")

    model$MIX <- as.data.frame(model$MIX)
    colnames(model$MIX) <- c("Iter", "Tree", "exp1", "tmin1", "tmax1", "exp2", "tmin2", "tmax2", "est", "kappa")

    model$mixNames <- c()
    for(i in 1:length(model$X)) {
      idx <- which(model$TreeStructs$exp == (i - 1)) # Converting indices from cpp to R

      if (length(idx) > 0) {
        model$TreeStructs$est[idx] <- model$TreeStructs$est[idx] * model$Yscale / model$X[[i]]$Xscale
      }

      for(j in i:length(model$X)) {
        if ((model$interaction > 1) || (j > i)) {
          model$mixNames <- c(model$mixNames, paste0(model$expNames[i], "-", model$expNames[j]))
        }

        idx <- which(model$MIX$exp1 == (i - 1) & model$MIX$exp2 == (j - 1))

        if (length(idx) > 0) {
          model$MIX$est[idx] <- model$MIX$est[idx] * model$Yscale / (model$X[[i]]$Xscale * model$X[[j]]$Xscale)
        }
      }
    }

    # Exposure names
    colnames(model$expProb) <- colnames(model$expCount) <- colnames(model$expInf) <- colnames(model$muExp) <- model$expNames
    if (model$interaction > 0) {
      colnames(model$mixInf) <- colnames(model$muMix) <- colnames(model$mixCount) <- model$mixNames
    }

    if (model$diagnostics) {
      model$treeAccept            <- as.data.frame(model$treeAccept)
      colnames(model$treeAccept)  <- c("tree", "step", "success", "exp", "term", "treeMhr", "mhr")
      model$treeAccept$step       <- factor(model$treeAccept$step,
                                            levels = 0:3,
                                            labels = c("Grow", "Prune", "Change", "SwitchExposure"))
      model$treeAccept$exp        <- factor(model$treeAccept$exp + 1,
                                            levels = 1:model$nExp,
                                            labels = model$expNames)
    }

  } else if (model$class %in% c("tdlm", "tdlnm", "monotone")) {
    # rescale DLM estimates
    model$TreeStructs           <- as.data.frame(model$TreeStructs)
    colnames(model$TreeStructs) <- c("Iter", "Tree", "xmin", "xmax", "tmin", "tmax", "est", "intcp")
    model$TreeStructs$est       <- model$TreeStructs$est * model$Yscale / model$Xscale
    model$TreeStructs$xmin      <- sapply(model$TreeStructs$xmin, function(i) {
                                              if (i == 0) {
                                                if (piecewise.linear) {min(model$X)}
                                                else {-Inf}
                                              } else {model$Xsplits[i]}
                                            }
                                          )
    model$TreeStructs$xmax      <- sapply(model$TreeStructs$xmax, function(i) {
                                              if (i == (length(model$Xsplits) + 1)) {
                                                if (piecewise.linear) {max(model$X)}
                                                else {Inf}
                                              } else {model$Xsplits[i]}
                                            }
                                          )
  } else {
    # Modifier output
    colnames(model$modProb) <- colnames(model$modCount) <- colnames(model$modInf) <- names(model$Mo)
    modNames    <- names(model$Mo)                     
    splitRules  <- strsplit(model$termRules, "&", TRUE)

    rule <- sapply(splitRules, function(str) {
              paste0(lapply(sort(str), function(rule) {
                # no rule
                if (length(rule) == 0) {
                  return("")
                # *** Continuous ***
                # >=
                } else if (length(spl <- strsplit(rule, ">=", TRUE)[[1]]) == 2) {
                  return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1], "']] >= ",
                                model$modSplitValRef[[as.numeric(spl[1]) + 1]][
                                  as.numeric(spl[2]) + 1]))
                # >
                } else if (length(spl <- strsplit(rule, "<", TRUE)[[1]]) == 2) {
                  return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1], "']] < ",
                                model$modSplitValRef[[as.numeric(spl[1]) + 1]][
                                  as.numeric(spl[2]) + 1]))
                # *** Categorical ***
                # in
                } else if (length(spl <- strsplit(rule, "[]", TRUE)[[1]]) == 2) {
                  inList <- paste0("c('", paste0(model$modSplitValRef[[
                    as.numeric(spl[1]) + 1]][
                      eval(parse(text = paste0("c(", spl[2], ")"))) + 1
                    ], collapse = "','"), "')")
                  return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1],
                                "']] %in% ", inList))
                # not in
                } else if (length(spl <- strsplit(rule, "][", TRUE)[[1]]) == 2) {
                  inList <- paste0("c('", paste0(model$modSplitValRef[[
                    as.numeric(spl[1]) + 1]][
                      eval(parse(text = paste0("c(", spl[2], ")"))) + 1
                    ], collapse = "','"), "')")
                  return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1],
                                "']] %notin% ", inList))
                } else {
                  return("")
                }
              }), collapse = " & ")})

    # Mixture interaction for hdlmm
    if (model$class == "hdlmm" & model$interaction > 0) {
      splitRulesMIX <- strsplit(model$termRuleMIX, "&", TRUE) 
      ruleMIX <- sapply(splitRulesMIX, function(str) {
                paste0(lapply(sort(str), function(rule) {
                  # no rule
                  if (length(rule) == 0) {
                    return("")
                  # *** Continuous ***
                  # >=
                  } else if (length(spl <- strsplit(rule, ">=", TRUE)[[1]]) == 2) {
                    return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1], "']] >= ",
                                  model$modSplitValRef[[as.numeric(spl[1]) + 1]][
                                    as.numeric(spl[2]) + 1]))
                  # >
                  } else if (length(spl <- strsplit(rule, "<", TRUE)[[1]]) == 2) {
                    return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1], "']] < ",
                                  model$modSplitValRef[[as.numeric(spl[1]) + 1]][
                                    as.numeric(spl[2]) + 1]))
                  # *** Categorical ***
                  # in
                  } else if (length(spl <- strsplit(rule, "[]", TRUE)[[1]]) == 2) {
                    inList <- paste0("c('", paste0(model$modSplitValRef[[
                      as.numeric(spl[1]) + 1]][
                        eval(parse(text = paste0("c(", spl[2], ")"))) + 1
                      ], collapse = "','"), "')")
                    return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1],
                                  "']] %in% ", inList))
                  # not in
                  } else if (length(spl <- strsplit(rule, "][", TRUE)[[1]]) == 2) {
                    inList <- paste0("c('", paste0(model$modSplitValRef[[
                      as.numeric(spl[1]) + 1]][
                        eval(parse(text = paste0("c(", spl[2], ")"))) + 1
                      ], collapse = "','"), "')")
                    return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1],
                                  "']] %notin% ", inList))
                  } else {
                    return("")
                  }
                }), collapse = " & ")})
    } else {
      ruleMIX = NA
    }
        
    model$modPairs <- sort(table(do.call(c, lapply(splitRules, function(r) {
      if (length(r) == 0) {
        return(NA)
      }

      m <- sapply(strsplit(r, ">=|<|\\[\\]|\\]\\["), function(i) as.numeric(i[1]))
      
      if (length(m) < 2) {
        return(modNames[m + 1])
      }
        
      c <- combn(length(m), 2)
      return(unique(sapply(1:ncol(c), function(i) paste0(modNames[sort(m[c[,i]]) + 1], collapse = "-"))))
    })))) / (model$nTrees * model$mcmcIter)


    # *** Combine the rules and the exposure data frames for HDLM, HDLMM ***
    if (model$class == "hdlmm") {
      # Combine DLM with rules with colnames
      model$TreeStructs           <- cbind.data.frame(rule, model$TreeStructs)
      colnames(model$TreeStructs) <- c("Rule", "Iter", "Tree", "Mod", "dlmPair", "dlmTerm", "exp", "tmin", "tmax", "est", "kappa")
      
      # Default of model$MIX is a vector of zeros
      if (model$interaction != 0) {
        model$MIX           <- as.data.frame(model$MIX)
        model$MIX           <- cbind.data.frame(ruleMIX, model$MIX)
        colnames(model$MIX) <- c("Rule", "Iter", "Tree", "Mod", "exp1", "tmin1", "tmax1", "exp2", "tmin2", "tmax2", "est")
      }

      model$mixNames <- c()

    } else {
      model$TreeStructs <- cbind.data.frame(rule, model$TreeStructs)
    }

    model$modSplitIdx <- NULL
    model$fullIdx     <- NULL


    # Redundant code
    if (model$class != "hdlmm") {
      colnames(model$TreeStructs) <- c("Rule", "Iter", "Tree", "modTerm", "dlnmTerm", "xmin", "xmax", "tmin", "tmax", "est")
      model$TreeStructs$est       <-  model$TreeStructs$est * model$Yscale / model$Xscale
      model$TreeStructs$xmin      <- sapply(model$TreeStructs$xmin, function(i) {
        if (i == 0) {-Inf}
        else {model$Xsplits[i]}
      })
      model$TreeStructs$xmax <- sapply(model$TreeStructs$xmax, function(i) {
        if (i == (length(model$Xsplits) + 1)) {Inf}
        else {model$Xsplits[i]}
      })
    } else { # hdlmm
      # Rescale: Iterate through modifiers and exposures to summarize TreeStruct and MIX
      for(i in 1:length(model$X)) {
        idx <- which(model$TreeStructs$exp == (i - 1)) # Converting indices from cpp to R
        if (length(idx) > 0)
          model$TreeStructs$est[idx] <- model$TreeStructs$est[idx] * model$Yscale / model$X[[i]]$Xscale

        for(j in i:length(model$X)) {
          if ((model$interaction > 1) || (j > i)){
            model$mixNames <- c(model$mixNames, paste0(model$expNames[i], "-", model$expNames[j]))
          }
            
          if (nrow(model$MIX) > 0) {
            idx <- which(model$MIX$exp1 == (i - 1) & model$MIX$exp2 == (j - 1))
            if (length(idx) > 0) {
              model$MIX$est[idx] <- model$MIX$est[idx] * model$Yscale / (model$X[[i]]$Xscale * model$X[[j]]$Xscale)
            }
          }
        }
      }
    }
  }

  # Remove model and exposure data unless stated otherwise
  model$data <- data
  if(!save.data){
    if (!(model$class %in% c("tdlmm", "hdlmm"))) {
      model$X     <- NULL
      model$Tcalc <- NULL
      model$Xcalc <- NULL
      model$Z     <- NULL
      model$Mo    <- NULL
    } else {
      for(i in 1:length(model$X)) {
        model$X[[i]]$X <- model$X[[i]]$Tcalc <- NULL
      }
      model$Z     <- NULL
      model$Z.zi  <- NULL
    }
  }
  
  # Change env to list
  model.out         <- lapply(names(model), function(i) model[[i]])
  names(model.out)  <- names(model)
  rm(model)
  gc()
  
  # return the call
  model.out$call <- match.call()
  
  class(model.out)  <- model.out$class
  return(model.out)
}