#' tdlmm
#' @title Treed Distributed Lag Mixture Models (Deprecated)
#' @description TDLMM is a method for estimating a Treed Distributed Lag
#' Mixture Model. It operates by building an ensemble of pairs of regression
#' trees. Each tree in a tree-pair partitions the time span of the exposure data
#' and estimates a piecewise constant distributed lag effect. The two trees
#' are then intersected to create an interaction surface for estimating the
#' interaction between two exposures. Exposures are selected for each tree
#' stochastically and each exposure or interaction has a unique shrinkage
#' variance component. This allows for exposure variable selection in addition
#' to the estimation of the distributed lag mixture model.
#'
#' @param formula object of class formula, a symbolic description of the fixed
#' effect model to be fitted, e.g. y ~ a + b
#' @param data data frame containing variables used in the formula
#' @param exposure.data named list containing equally sized numerical matrices
#' of exposure data with same, having same length as data
#' @param n.trees integer for number of trees in ensemble
#' @param n.burn integer for length of burn-in
#' @param n.iter integer for number of iterations to run model after burn-in
#' @param n.thin integer thinning factor, i.e. keep every tenth iteration
#' @param family 'gaussian' for continuous response, 'logit' for binomial
#' response with logit link, or 'zinb' for zero-inflated negative binomial with logit link
#' @param binomial.size integer type scalar (if all equal, default = 1) or
#' vector defining binomial size for 'logit' family
#' @param formula.zi object of class formula, a symbolic description of the ZI
#' model to be fitted, e.g. y ~ a + b. This only applies to ZINB where covariates for
#' ZI model is different from NB model. This is same as the main formula by default
#' @param keep_XZ FALSE (default) or TRUE: keep the model scale exposure and covariate data
#' @param mixture.interactions 'noself' (default) which estimates interactions
#' only between two different exposures, 'all' which also allows
#' interactions within the same exposure, or 'none' which eliminates all
#' interactions and estimates only main effects of each exposure
#' @param tree.params numerical vector of alpha and beta hyperparameters
#' controlling tree depth (see Bayesian CART, 1998), default: alpha = 0.95,
#' beta = 2
#' @param step.prob numerical vector for probability of 1) grow/prune,
#' 2) change, 3) switch exposure, defaults to (0.25, 0.25, 0.25) or equal
#' probability of each step for tree updates
#' @param mix.prior positive scalar hyperparameter for sparsity of exposures
#' @param shrinkage character "all" (default), "trees", "exposures", "none",
#' turns on horseshoe-like shrinkage priors for different parts of model
#' @param subset integer vector to analyze only a subset of data and exposures
#' @param verbose TRUE (default) or FALSE: print output
#' @param diagnostics TRUE or FALSE (default) keep model diagnostic such as
#' terminal nodes, acceptance details, etc.
#' @param initial.params initial parameters for fixed effects model, FALSE = none (default), "glm" = generate using GLM, or user defined, length must equal number of parameters in fixed effects model
#' @param ... NA
#'
#' @details Model is recommended to be run for at minimum 5000 burn-in
#' iterations followed by 15000 sampling iterations with a thinning factor of 5.
#' Convergence can be checked by re-running the model and validating consistency
#' of results.
#'
#' @returns object of class 'tdlmm'
#'
tdlmm <- function(formula,
                  data,
                  exposure.data,
                  n.trees = 20,
                  n.burn = 2000,
                  n.iter = 5000,
                  n.thin = 5,
                  family = "gaussian",
                  binomial.size = 1,
                  formula.zi = NULL, 
                  keep_XZ = FALSE,
                  mixture.interactions = "noself",
                  tree.params = c(.95, 2),
                  step.prob = c(.25, .25, .25),
                  mix.prior = 1,
                  shrinkage = "exposures",
                  subset = NULL,
                  verbose = TRUE,
                  diagnostics = FALSE,
                  initial.params = NULL,
                  ...)
{
  message("Execution halt: 'tdlmm.R' has been deprecated and will be removed in the future. Please use 'dlmtree.fit' with arguments: dlm.type = `linear` and mixture = TRUE, instead. For an example code, use ?dlmtree.fit.")

  # model <- list()                 
  # options(stringsAsFactors = F)

  # # ---- Check user inputs ----
  # # number of iterations [Control for MCMC iteration vs thinning]
  # if (n.iter < n.thin * 10)
  #   stop("after thinning you will be left with less than 10 MCMC samples,
  #        increase the number of iterations!")
  
  # # data for formula [data frame check]
  # if (!is.data.frame(data))
  #   stop("`data` must be a data.frame")

  # # exposure data
  # if (!is.list(exposure.data))
  #   stop("`exposure.data` must be a list of named exposures")

  # model$nExp <- length(exposure.data)  
  # model$expNames <- names(exposure.data)

  # # [Exposure name check]
  # if (is.null(model$expNames) || length(unique(model$expNames)) != model$nExp || any(model$expNames == ""))
  #   stop("`exposure.data` must be a named list with unique, non-empty names")

  # model$pExp <- ncol(exposure.data[[1]])

  # # Sanity check for each exposure 
  # for (i in 1:length(exposure.data)) {
  #   # [Not numeric values for exposure data]
  #   if (!is.numeric(exposure.data[[i]]))
  #     stop("each exposure in list `exposure.data` must be a numeric matrix")
  #   # [Number of observations do not match]
  #   if (nrow(data) != nrow(exposure.data[[i]]))
  #     stop("`data` and `exposure.data` must have same number of rows")
  #   # [Exposure total time lag does not match with the others]
  #   if (ncol(exposure.data[[i]]) != model$pExp)
  #     stop("each exposure in `exposure.data` must have the same
  #          number of time observations")
  #   # [NA/Missing values in exposure data]
  #   if (any(is.na(exposure.data[[i]]))) # any(): Is there at least one which meets the condition?
  #     stop("missing values in exposure data")
  # }

  # # [iteration control] Check if all MCMC parameters are positive integers
  # if (all(sapply(list(n.trees, n.burn, n.iter, n.thin),
  #                function(i) is.integer(i) & i > 0)))
  #   stop("n.* must be integer and > 0")

  # # [response type] Model much be Gaussian / Logistic / Zero-inflated negative binomial
  # if (!(family %in% c("gaussian", "logit", "zinb"))) 
  #   stop("`family` must be one of `gaussian`, `logit`, or 'zinb'")

  # # [binomial size] Only for binomial 
  # model$binomial <- 0                     
  # model$binomialSize <- rep(0, nrow(data))

  # # If binomial is called, update the parameters(binomial, binomialSize) above.
  # if (family == "logit") {
  #   # If binary, update binomialSize to a vector of ones.
  #   if (length(binomial.size) == 1)                       
  #     binomial.size <- rep(binomial.size, nrow(data)) 

  #   # binomialSize vector must be the same length as the sample size.
  #   if (length(binomial.size) != nrow(data))
  #     stop("`binomial.size` must be positive integer and same length as data")

  #   # Update the parameters
  #   model$binomialSize <- force(binomial.size)
  #   model$binomial <- 1
  # }

  # # ZINB flag
  # model$zinb <- 0
  # if (family == "zinb") {
  #   model$zinb <- 1
  #   model$sigma2 <- 1
  # }

  # # [Mixture interactions]
  # # [Correct options for interaction parameter]
  # if (!(mixture.interactions %in% c("noself", "all", "none")))
  #   stop("`mixture.interactions must be one of `noself`, `all`, `none`")

  # if (mixture.interactions %in% c("marginal", "none")) {
  #   model$interaction <- 0     
  #   model$nMix <- 0   
  # } else if (mixture.interactions == "noself") { 
  #   model$interaction <- 1  
  #   model$nMix <- model$nExp * (model$nExp - 1) / 2 
  # } else { 
  #   model$interaction <- 2   
  #   model$nMix <- model$nExp * (model$nExp + 1) / 2 
  # }

  # # [tree parameters]
  # # Tree param default: alpha, beta = c(.95, 2)
  # if (length(tree.params) != 2) {
  #   stop("tree.params must have length 2")
  # } else {
  #   # [Sanity check on alpha & beta range]
  #   if (tree.params[1] > 1 | tree.params[1] < 0)
  #     stop("tree.params[1] must be between 0 and 1")
  #   if (tree.params[2] < 0)
  #     stop("tree.params[2] must be greater than 0")
  # }

  # # step probabilities
  # # 1) grow/prune, 2) change, 3) switch exposure, defaults to (0.25, 0.25, 0.25)
  # if (any(step.prob < 0) || any(step.prob > 1) || length(step.prob) != 3)
  #   stop("must specify three step.prob components between zero and one")


  # # ---- [Model control arguments] ----
  # model$nTrees <- force(n.trees) 
  # model$nBurn <- force(n.burn)
  # model$nIter <- force(n.iter) 
  # model$nThin <- force(n.thin)
  # model$mcmcIter <- force(floor(n.iter / n.thin))
  # model$family <- force(family) 
  # model$verbose <- force(verbose) 
  # model$diagnostics <- force(diagnostics) 
  # model$treePriorTDLM <- force(tree.params) 
  
  # # Step probabilities
  # model$stepProb <- force(c(step.prob[1], step.prob[1], step.prob[2], step.prob[3])) 
  # model$stepProbTDLM <- force(model$stepProb / sum(model$stepProb))
  # model$mixPrior <- ifelse(mix.prior == "auto", -1, mix.prior)
  # model$shrinkage <- ifelse(shrinkage == "all", 3,
  #                           ifelse(shrinkage == "trees", 2,
  #                                  ifelse(shrinkage == "exposures", 1, 0)))

  # if (model$verbose)
  #   cat("Preparing data...\n")

  # # ---- [Create data subset] ----
  # if (!is.null(subset)) {
  #   if (length(subset) > 1 & is.integer(subset) & 
  #       all(subset > 0) & all(subset <= nrow(data))) {
  #     data <- data[subset,]
  #     exposure.data <- lapply(exposure.data, function(i) i[subset,]) 
  #     if (model$family == "logit")
  #       model$binomialSize <- model$binomialSize[subset]
  #   } else {
  #     stop("invalid subset, must be integers within range of data length")
  #   }
  # }


  # # ---- [Setup control and response variables] ----
  # if (is.null(formula.zi)) {
  #   formula.zi <- formula
  # }

  # model$formula <- force(as.formula(formula))     
  # model$formula.zi <- force(as.formula(formula.zi)) 
  # tf <- terms.formula(model$formula, data = data) 
  # tf.zi <- terms.formula(model$formula.zi, data = data) 

  # # Sanity check for response
  # if (!attr(tf, "response") & !attr(tf.zi, "response"))
  #   stop("no valid response in formula")

  # # Save if intercept is included in the model
  # model$intercept <- force(ifelse(attr(tf, "intercept"), TRUE, FALSE)) 
  # model$intercept.zi <- force(ifelse(attr(tf.zi, "intercept"), TRUE, FALSE)) 

  # # Sanity check for variables and the names
  # if (length(which(attr(tf, "term.labels") %in% colnames(data))) == 0 & 
  #     length(which(attr(tf.zi, "term.labels") %in% colnames(data))) == 0 & 
  #     !model$intercept & !model$intercept.zi) 
  #   stop("no valid variables in formula")

  # # ---- [Exposure splits] ----
  # model$splitProb <- as.double(c())
  # model$timeProb <- force(rep(1 / (model$pExp - 1), model$pExp - 1))
  
  # # Exposure data organize
  # model$X <- list() 
  # for (i in 1:model$nExp) {
  #   model$X[[i]] <- force(list(Xscale = sd(exposure.data[[i]]),   
  #                              X = exposure.data[[i]]))           
  #   model$X[[i]]$X <- force(model$X[[i]]$X / model$X[[i]]$Xscale) 
  #   model$X[[i]]$Xrange <- force(range(model$X[[i]]$X))           
  #   model$X[[i]]$Xquant <- force(quantile(model$X[[i]]$X, 0:100/100) *  model$X[[i]]$Xscale)
  #   model$X[[i]]$intX <- force(mean(model$X[[i]]$X))              
  #   model$X[[i]]$Tcalc <-
  #     force(sapply(1:ncol(model$X[[i]]$X), function(j) { 
  #       rowSums(model$X[[i]]$X[, 1:j, drop = F]) }))       
  # }
  # names(model$X) <- model$expNames                                
  # model$expProb <- force(rep(1/length(model$X), length(model$X))) 


  # # ---- [Scale data and setup exposures] ----
  # data <- droplevels(data)                                     
  # mf <- model.frame(model$formula, data = data,                
  #                   drop.unused.levels = TRUE,
  #                   na.action = NULL)
  # mf.zi <- model.frame(model$formula.zi, data = data,               
  #                     drop.unused.levels = TRUE,
  #                     na.action = NULL)
  # if (any(is.na(mf)) & any(is.na(mf.zi)))
  #   stop("missing values in model data, use `complete.cases()` to subset data")


  # # Dan's ---
  # # model$Y <- force(model.response(mf))
  # # model$Z <- force(model.matrix(model$formula, data = mf))
  # # # QR <- qr(crossprod(model$Z))
  # # # model$Z <- model$Z[,sort(QR$pivot[seq_len(QR$rank)])]
  # # model$Znames <- colnames(model$Z)#[sort(QR$pivot[seq_len(QR$rank)])]
  # # model$droppedCovar <- c()#colnames(model$Z)[QR$pivot[-seq_len(QR$rank)]]
  # # # model$Z <- matrix(model$Z[,sort(QR$pivot[seq_len(QR$rank)])], nrow(model$Z), QR$rank)
  # # # model$Z <- force(scaleModelMatrix(model$Z))
  # # # rm(QR)
  
  # # Organize response variable & fixed effect variable
  # # Main covariates
  # model$Y <- force(model.response(mf))                        
  # model$Z <- force(model.matrix(model$formula, data = mf))  
  # QR <- qr(crossprod(model$Z))                              
  # model$Z <- model$Z[, sort(QR$pivot[seq_len(QR$rank)])]      
  # model$Znames <- colnames(model$Z)[sort(QR$pivot[seq_len(QR$rank)])]
  # model$droppedCovar <- colnames(model$Z)[QR$pivot[-seq_len(QR$rank)]]
  # model$Z <- force(scaleModelMatrix(model$Z))
  # rm(QR)

  # # covariates for ZI model
  # model$Z.zi <- force(model.matrix(model$formula.zi, data = mf.zi))  
  # QR.zi <- qr(crossprod(model$Z.zi))                                 
  # model$Z.zi <- model$Z.zi[,sort(QR.zi$pivot[seq_len(QR.zi$rank)])]           
  # model$droppedCovar.zi <- colnames(model$Z.zi)[QR.zi$pivot[-seq_len(QR.zi$rank)]]
  # model$Z.zi <- force(scaleModelMatrix(model$Z.zi))
  # rm(QR.zi)

  # # Organize for different models
  # # Gaussian
  # if (model$family == "gaussian") {
  #   model$Ymean <- sum(range(model$Y)) / 2
  #   model$Yscale <- diff(range(model$Y - model$Ymean))
  #   model$Y <- force((model$Y - model$Ymean) / model$Yscale)
  # } else { # For both logistic and ZINB,
  #   model$Yscale <- 1
  #   model$Ymean <- 0
  #   model$Y <- force(scale(model$Y, center = 0, scale = 1))
  # }

  # # Store the processed values
  # model$Y <- force(c(model$Y))

  # model$Zscale <- attr(model$Z, "scaled:scale")
  # model$Zmean <- attr(model$Z, "scaled:center")
  # model$Znames <- colnames(model$Z)
  # model$Z <- force(matrix(model$Z, nrow(model$Z), ncol(model$Z)))

  # model$Zscale.zi <- attr(model$Z.zi, "scaled:scale")
  # model$Zmean.zi <- attr(model$Z.zi, "scaled:center")
  # model$Znames.zi <- colnames(model$Z.zi)
  # model$Z.zi <- force(matrix(model$Z.zi, nrow(model$Z.zi), ncol(model$Z.zi)))

  
  
  
  
  # model$initParams <- rep(0, ncol(model$Z))
  # if (!is.null(initial.params)) {
  #   names(model$initParams) <- colnames(model$Z)
  #   if (sum(names(initial.params) %in% colnames(model$Z)) > 0) {
  #     na <- names(initial.params[ # get matching names 
  #       which(names(initial.params) %in% colnames(model$Z))])
  #     model$initParams[na] <- initial.params[na]
  #   }
  # }

  # # ---- Run model ----
  # out <- tdlmm_Cpp(model)

  # if (verbose)
  #   cat("\nCompiling results...\n")

  # for (n in names(out))
  #   model[[n]] <- out[[n]]



  # # ---- Prepare output ----
  # model$Y <- model$Y * model$Yscale + model$Ymean
  # model$fhat <- model$fhat * model$Yscale
  # model$sigma2 <- model$sigma2 * (model$Yscale^2)
  
  # # rescale fixed effect estimates
  # # Gaussian & Logistic fixed effect
  # model$gamma <- sapply(1:ncol(model$gamma), function(i) {
  # model$gamma[,i] * model$Yscale / model$Zscale[i] })

  # # ZINB fixed effect (binary)
  # model$b1 <- sapply(1:ncol(model$b1), function(i) {
  # model$b1[,i] * model$Yscale / model$Zscale.zi[i] })

  # # ZINB fixed effect (NegBin)
  # model$b2 <- sapply(1:ncol(model$b2), function(i) {
  # model$b2[,i] * model$Yscale / model$Zscale[i] })

  # # If intercept:
  # if (model$intercept) {
  #   model$gamma[,1] <- model$gamma[,1] + model$Ymean
  #   if (ncol(model$Z) > 1)
  #     model$gamma[,1] <- model$gamma[,1] - model$gamma[,-1,drop=FALSE] %*% model$Zmean[-1]

  #   model$b2[,1] <- model$b2[,1] + model$Ymean
  #   if (ncol(model$Z) > 1)
  #     model$b2[,1] <- model$b2[,1] - model$b2[,-1,drop=FALSE] %*% model$Zmean[-1]
  # }

  # if (model$intercept.zi) {
  #   model$b1[,1] <- model$b1[,1] + model$Ymean
  #   if (ncol(model$Z.zi) > 1)
  #     model$b1[,1] <- model$b1[,1] - model$b1[,-1,drop=FALSE] %*% model$Zmean.zi[-1]
  # }

  # colnames(model$gamma) <- model$Znames
  # colnames(model$b1) <- model$Znames.zi
  # colnames(model$b2) <- model$Znames

  # # rescale DLM and Mixture estimates
  # model$DLM <- as.data.frame(model$DLM)
  # colnames(model$DLM) <- c("Iter", "Tree", "TreePair", "exp", "tmin", "tmax",
  #                          "est", "kappa")
  # model$MIX <- as.data.frame(model$MIX)
  # colnames(model$MIX) <- c("Iter", "Tree", "exp1", "tmin1", "tmax1",
  #                          "exp2", "tmin2", "tmax2", "est", "kappa")
  # model$mixNames <- c()

  # for (i in 1:length(model$X)) {
  #   idx <- which(model$DLM$exp == (i - 1)) # Converting indices from cpp to R
  #   if (length(idx) > 0)
  #     model$DLM$est[idx] <- model$DLM$est[idx] * model$Yscale / model$X[[i]]$Xscale

  #   for (j in i:length(model$X)) {

  #     if ((model$interaction > 1) || (j > i))
  #       model$mixNames <- c(model$mixNames,
  #                           paste0(model$expNames[i], "-", model$expNames[j]))

  #     idx <- which(model$MIX$exp1 == (i - 1) & model$MIX$exp2 == (j - 1))
  #     if (length(idx) > 0) {
  #       model$MIX$est[idx] <- model$MIX$est[idx] * model$Yscale / (model$X[[i]]$Xscale * model$X[[j]]$Xscale)
  #     }
  #   }
  # }

  # # Exposure names
  # colnames(model$expProb) <- colnames(model$expCount) <-
  #   colnames(model$expInf) <- colnames(model$muExp) <- model$expNames
  # if (model$interaction > 0) {
  #   colnames(model$mixInf) <- colnames(model$muMix) <-
  #     colnames(model$mixCount) <- model$mixNames
  # }

  # if (model$diagnostics) {
  #   model$treeAccept <- as.data.frame(model$treeAccept)
  #   colnames(model$treeAccept) <- c("tree", "step", "success", "exp", "term", "treeMhr", "mhr")
  #   model$treeAccept$step <- factor(model$treeAccept$step,
  #                                   levels = 0:3,
  #                                   labels = c("Grow", "Prune", "Change",
  #                                              "SwitchExposure"))
  #   model$treeAccept$exp <- factor(model$treeAccept$exp + 1,
  #                                  levels = 1:model$nExp,
  #                                  labels = model$expNames)
  # }

  # # Remove model and exposure data
  # if (!keep_XZ) {
  #   for (i in 1:length(model$X))
  #     model$X[[i]]$X <- model$X[[i]]$Tcalc <- NULL
  #   model$Z <- NULL
  #   model$Z.zi <- NULL
  # }


  # # Change env to list
  # model.out <- lapply(names(model), function(i) model[[i]])
  # names(model.out) <- names(model)
  # rm(model)
  # gc()

  model.out <- NULL
  #class(model.out) <- "tdlmm"
  return(model.out)
}
