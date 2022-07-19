#' tdlmm
#' @title Treed Distributed Lag Mixture Models
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
#' response with logit link, or 'zinb' for zero-inflated negative binomial model
#' @param binomial.size integer type scalar (if all equal, default = 1) or
#' vector defining binomial size for 'logit' family
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
#' @param verbose true (default) or false: print output
#' @param diagnostics true or false (default) keep model diagnostic such as
#' terminal nodes, acceptance details, etc.
#' @param ... NA
#'
#' @details Model is recommended to be run for at minimum 5000 burn-in
#' iterations followed by 15000 sampling iterations with a thinning factor of 5.
#' Convergence can be checked by re-running the model and validating consistency
#' of results.
#'
#' @return object of class 'tdlmm'
#' @export
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
                  mixture.interactions = "noself",
                  tree.params = c(.95, 2),
                  step.prob = c(.25, .25, .25),
                  mix.prior = 1,
                  shrinkage = "all",
                  subset = NULL,
                  verbose = TRUE,
                  diagnostics = FALSE,
                  spNodes = NULL,
                  areaW = NULL,
                  areaA = NULL,
                  #geoid = NULL,
                  ...)
{
  model <- list()                 # Create a list to send to tdlmm_Cpp.cpp 
  options(stringsAsFactors = F)   # String is not a factor

  # ---- Check user inputs ----
  # number of iterations [Control for MCMC iteration vs thinning]
  if (n.iter < n.thin * 10)
    stop("after thinning you will be left with less than 10 MCMC samples,
         increase the number of iterations!")
  
  # data for formula [data frame check]
  if (!is.data.frame(data))
    stop("`data` must be a data.frame")

  # exposure data
  if (!is.list(exposure.data))
    stop("`exposure.data` must be a list of named exposures")

  model$nExp <- length(exposure.data)     # Number of exposure
  model$expNames <- names(exposure.data)  # Names of exposures 

  # [Exposure name check] exposure names cannot be empty & # of exposure names must match with # of exposures
  if (is.null(model$expNames) || length(unique(model$expNames)) != model$nExp || any(model$expNames == ""))
    stop("`exposure.data` must be a named list with unique, non-empty names")

  model$pExp <- ncol(exposure.data[[1]])  # Total time lag (37)

  # Sanity check for each exposure 
  for (i in 1:length(exposure.data)) {    # length(exposure.data) = Number of exposures (nExp)
    # [Not numeric values for exposure data]
    if (!is.numeric(exposure.data[[i]]))
      stop("each exposure in list `exposure.data` must be a numeric matrix")
    # [Number of observations do not match]
    if (nrow(data) != nrow(exposure.data[[i]]))
      stop("`data` and `exposure.data` must have same number of rows")
    # [Exposure total time lag does not match with the others]
    if (ncol(exposure.data[[i]]) != model$pExp)
      stop("each exposure in `exposure.data` must have the same
           number of time observations")
    # [NA/Missing values in exposure data]
    if (any(is.na(exposure.data[[i]]))) # any(): Is there at least one which meets the condition?
      stop("missing values in exposure data")
  }

  # [iteration control] Check if all MCMC parameters are positive integers
  if (all(sapply(list(n.trees, n.burn, n.iter, n.thin), # all() = opposite of any()
                 function(i) is.integer(i) & i > 0)))
    stop("n.* must be integer and > 0")

  # [response type] Model much be Gaussian / Logistic / Zero-inflated negative binomial
  if (!(family %in% c("gaussian", "logit", "zinb"))) 
    stop("`family` must be one of `gaussian`, `logit`, or 'zinb'")

  # [binomial size] Only for binomial 
  # If not binomial, 
  model$binomial <- 0                       # Set binomial flag to 0.
  model$binomialSize <- rep(0, nrow(data))  # binomialSize is a vector of zero 

  # If binomial is called, update the parameters(binomial, binomialSize) above.
  if (family == "logit") {
    # If binary, update binomialSize to a vector of ones.
    if (length(binomial.size) == 1)                       
      binomial.size <- rep(binomial.size, nrow(data)) 

    # binomialSize vector must be the same length as the sample size.
    if (length(binomial.size) != nrow(data))
      stop("`binomial.size` must be positive integer and same length as data")

    # Update the parameters
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

  # # geoid check -> UPDATE WHEN WE HAVE AN ACTUAL DATASET WITH GEOID TYPE DATA
  # if(!is.null(geoid)){
  #   if(nrow(data) != length(geoid)){ 
  #     stop("every observation must have a geoid")
  #   }

  #   if(sum(nchar(geoid)) != nrow(data)*11){
  #     stop("geoid must have 11 letters")
  #   }
  # }

  # [Mixture interactions]
  # [Correct options for interaction parameter]
  if (!(mixture.interactions %in% c("noself", "all", "none")))
    stop("`mixture.interactions must be one of `noself`, `all`, `none`")

  if (mixture.interactions %in% c("marginal", "none")) { # TDLMMadd
    model$interaction <- 0      # No interaction
    model$nMix <- 0             # No interaction terms
  } else if (mixture.interactions == "noself") { # TDLMMns
    model$interaction <- 1      # All interaction except oneself
    model$nMix <- model$nExp * (model$nExp - 1) / 2 # (M choose 2) = M(M-1)/2
  } else {  # TDLMMall
    model$interaction <- 2      # All interaction including oneself
    model$nMix <- model$nExp * (model$nExp + 1) / 2 # (M choose 2) + M(exposure itself) = M(M+1)/2
  }

  # [tree parameters]
  # Tree param default: alpha, beta = c(.95, 2)
  if (length(tree.params) != 2) {
    stop("tree.params must have length 2")
  } else {
    # [Sanity check on alpha & beta range]
    if (tree.params[1] > 1 | tree.params[1] < 0)
      stop("tree.params[1] must be between 0 and 1")
    if (tree.params[2] < 0)
      stop("tree.params[2] must be greater than 0")
  }

  # step probabilities
  # 1) grow/prune, 2) change, 3) switch exposure, defaults to (0.25, 0.25, 0.25)
  if (any(step.prob < 0) || any(step.prob > 1) || length(step.prob) != 3)
    stop("must specify three step.prob components between zero and one")


  # ---- [Model control arguments] ----
  model$nTrees <- force(n.trees) # Number ot trees
  model$nBurn <- force(n.burn)   # Number of burn-in
  model$nIter <- force(n.iter)   # Number of iteration
  model$nThin <- force(n.thin)   # Number of thinning
  model$mcmcIter <- force(floor(n.iter / n.thin)) # rounding MCMC iteration / thinning
  model$family <- force(family)  # Gaussian / Logit / ZINB
  model$verbose <- force(verbose) # Boolean for print output
  model$diagnostics <- force(diagnostics) # keep model diagnostic
  model$treePrior <- force(tree.params)   # alpha, beta
  model$stepProb <- force(c(step.prob[1], step.prob[1], step.prob[2], step.prob[3])) # Step probabilities
  model$stepProb <- force(model$stepProb / sum(model$stepProb))
  model$mixPrior <- mix.prior # Positive scalar hyperparameter for sparsity of exposures (Linero)
  model$shrinkage <- ifelse(shrinkage == "all", 3,
                            ifelse(shrinkage == "trees", 2,
                                   ifelse(shrinkage == "exposures", 1, 0))) # Assign shrinkage to 3, 2, 1 (horseshoe)

  if (model$verbose)
    cat("Preparing data...\n")

  # ---- [Create data subset] ----
  if (!is.null(subset)) {
    if (length(subset) > 1 & is.integer(subset) & # If subset index vector is not null, all the indices are proper,
        all(subset > 0) & all(subset <= nrow(data))) {
      data <- data[subset,] # Subset the data
      exposure.data <- lapply(exposure.data, function(i) i[subset,]) # Subset the rows of exposure data
      if (model$family == "logit")
        model$binomialSize <- model$binomialSize[subset]
    } else {
      stop("invalid subset, must be integers within range of data length")
    }
  }


  # ---- [Setup control and response variables] ----
  model$formula <- force(as.formula(formula))     # Save the string "y ~ ." as a formula type
  tf <- terms.formula(model$formula, data = data) # terms.formula with data -> returns model.matrix & attributes
  
  # Sanity check for response
  if (!attr(tf, "response"))
    stop("no valid response in formula")

  # Save if intercept is included in the model
  model$intercept <- force(ifelse(attr(tf, "intercept"), TRUE, FALSE)) 

  # Sanity check for variables and the names
  if (length(which(attr(tf, "term.labels") %in% colnames(data))) == 0) 
    stop("no valid variables in formula")

  # ---- [Exposure splits] ----
  model$splitProb <- as.double(c())
  model$timeProb <- force(rep(1 / (model$pExp - 1), model$pExp - 1)) # Uniform probability 1/36 of 36 - time lags
  
  # Exposure data organize
  model$X <- list() 
  for (i in 1:model$nExp) { # For each exposure,
    model$X[[i]] <- force(list(Xscale = sd(exposure.data[[i]]),     # Store standard deviation of exposure data
                               X = exposure.data[[i]]))             # Store exposure data itself
    model$X[[i]]$X <- force(model$X[[i]]$X / model$X[[i]]$Xscale)   # Scale the stored exposure data
    model$X[[i]]$Xrange <- force(range(model$X[[i]]$X))             # Store range of exposure data
    model$X[[i]]$Xquant <- force(quantile(model$X[[i]]$X, 0:100/100) *  model$X[[i]]$Xscale) # Store unscaled quantile
    model$X[[i]]$intX <- force(mean(model$X[[i]]$X))                # Store the mean of scaled exposure data
    # Tcalc is a matrix (n x T) where column j, (j = 1,...,T) represents a rowsum of exposure data from t = 1 to t = j
    model$X[[i]]$Tcalc <-
      force(sapply(1:ncol(model$X[[i]]$X), function(j) { 
        rowSums(model$X[[i]]$X[, 1:j, drop = F]) }))       
  }
  names(model$X) <- model$expNames                                  # Save exposure names
  model$expProb <- force(rep(1/length(model$X), length(model$X)))   # Uniform probability for exposure selection


  # ---- [Scale data and setup exposures] ----
  data <- droplevels(data)                                        # Drop all levels of factors in the data.
  mf <- model.frame(model$formula, data = data,                   # Extract Model information with formula and data
                    drop.unused.levels = TRUE,
                    na.action = NULL)
  if (any(is.na(mf)))
    stop("missing values in model data, use `complete.cases()` to subset data")
  
  # Organize response variable & fixed effect variable
  model$Y <- force(model.response(mf))                            # Response Y
  model$Z <- force(model.matrix(model$formula, data = mf))        # Fixed effect matrix, Z (c1 ~ c5, b1 ~ b5)
  QR <- qr(crossprod(model$Z))                                    # t(Z)%*%Z and its QR decomposition
  model$Z <- model$Z[,sort(QR$pivot[seq_len(QR$rank)])]           
  model$droppedCovar <- colnames(model$Z)[QR$pivot[-seq_len(QR$rank)]]
  model$Z <- force(scaleModelMatrix(model$Z))
  rm(QR)

  # Organize for different models
  # Gaussian
  if (model$family == "gaussian") {
    model$Ymean <- sum(range(model$Y)) / 2
    model$Yscale <- diff(range(model$Y - model$Ymean))
    model$Y <- force((model$Y - model$Ymean) / model$Yscale)
  } else { # For both logistic and ZINB,
    model$Yscale <- 1
    model$Ymean <- 0
    model$Y <- force(scale(model$Y, center = 0, scale = 1))
  }

  # Store the processed values
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

  # FOR REAL DATASET
  #if(!is.null(geoid)){
    # if (model$verbose){
    #   cat("Gathering spatial data...\n")
    # }

    # stateCode <- substr(geoid[1], 1, 2)          # Could be a vector but usually one number
    # countiesCode <- unique(substr(geoid, 3, 5))  # A vector

    # # With geoid, create a SpatialPolygonsDataFrame (SPDF) object
    # ctspatial <- tigris::tracts(state = stateCode, county = countiesCode, refresh = T, class = "sp") 

    # # geoid_df & sp_data
    # geoid_df = data.frame("geoid" = geoid)
    # geoid_df$state = substr(geoid_df$geoid, 1, 2)  
    # geoid_df$county = substr(geoid_df$geoid, 3, 5)  
    # geoid_df$tract = substr(geoid_df$geoid, 6, 11)  

    # sp_data = cbind(data, geoid_df)

    # # Merge the spatial information with the cityfile
    # spatial_data <- sp::merge(ctspatial, sp_data, by.y = "geoid", by.x = "GEOID", all.x = F, all.y = T, duplicateGeoms = TRUE)

    # # Get an adjacency matrix
    # adj <- poly2nb(spatial_data)

    # #adj <- poly2nb(spatial_data, row.names = spatial_data$GEOID)
    # coord_df = coordinates(spatial_data)
    # row.names(coord_df) = spatial_data$GEOID

    # W <- nb2mat(adj, style = "B") # Adjacency matrix
    # D <- diag(rowSums(W))         # Diagonal matrix

    # # Compute variance-covariance matrix for phi (CAR model)
    # alpha = 0.5 # alpha = 0 implies no spatial effect. Note that alpha = 1 makes Q non-invertible
    # Q <- chol2inv(chol(D - alpha*W))
    # model$Q = Q
    # model$spatial = TRUE
  #}

  # ---- Run model ----
  # model.list <- lapply(ls(envir = model), function(i) model[[i]])
  # names(model.list) <- ls(envir = model)
  out <- tdlmm_Cpp(model)

  # if (model$family == "gaussian")
  #   out <- tdlmmGaussian(model)
  # else
  #   out <- tdlmmBinomial(model)

  # rm(model.list)
  if (verbose)
    cat("\nCompiling results...\n")

  for (n in names(out))
    model[[n]] <- out[[n]]


  # ---- Prepare output ----
  model$Y <- model$Y * model$Yscale + model$Ymean
  model$fhat <- model$fhat * model$Yscale
  model$sigma2 <- model$sigma2 * (model$Yscale ^ 2)
  
  # rescale fixed effect estimates
  #if (model$intercept) {
  #  model$gamma[,-1] <- sapply(2:ncol(model$gamma), function(i) { # Logistic
  #    model$gamma[,i] * model$Yscale / model$Zscale[i]})
  #  model$b1[,-1] <- sapply(2:ncol(model$b1), function(i) { # ZINB binary
  #    model$b1[,i] * model$Yscale / model$Zscale[i]})
  #  model$b2[,-1] <- sapply(2:ncol(model$b2), function(i) { # ZINB count
  #    model$b2[,i] * model$Yscale / model$Zscale[i]})
  #} else {
  #  model$gamma <- sapply(1:ncol(model$gamma), function(i) { # Logistic
  #    model$gamma[,i] * model$Yscale / model$Zscale[i]})
  #  model$b1 <- sapply(1:ncol(model$b1), function(i) { # ZINB binary
  #    model$b1[,i] * model$Yscale / model$Zscale[i]})
  #  model$b2 <- sapply(1:ncol(model$b2), function(i) { # ZINB count
  #    model$b2[,i] * model$Yscale / model$Zscale[i]})
  #}

  # Gaussian & Logistic fixed effect
  model$gamma <- sapply(1:ncol(model$gamma), function(i) {
  model$gamma[,i] * model$Yscale / model$Zscale[i] })

  # ZINB fixed effect (binary)
  model$b1 <- sapply(1:ncol(model$b1), function(i) {
  model$b1[,i] * model$Yscale / model$Zscale[i] })

  # ZINB fixed effect (NegBin)
  model$b2 <- sapply(1:ncol(model$b2), function(i) {
  model$b2[,i] * model$Yscale / model$Zscale[i] })

  # ZINB + Spatial random effect
  model$spPhi <- sapply(1:ncol(model$spPhi), function(i) {
  model$spPhi[,i] * model$Yscale / model$Ascale[i] })

  if (model$intercept) {
    model$gamma[,1] <- model$gamma[,1] + model$Ymean
    if (ncol(model$Z) > 1)
      model$gamma[,1] <- model$gamma[,1] - model$gamma[,-1] %*% model$Zmean[-1]

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

  # rescale DLM and Mixture estimates
  model$DLM <- as.data.frame(model$DLM)
  colnames(model$DLM) <- c("Iter", "Tree", "exp", "tmin", "tmax",
                           "est", "kappa")
  model$MIX <- as.data.frame(model$MIX)
  colnames(model$MIX) <- c("Iter", "Tree", "exp1", "tmin1", "tmax1",
                           "exp2", "tmin2", "tmax2", "est", "kappa")
  model$mixNames <- c()

  for (i in 1:length(model$X)) {

    idx <- which(model$DLM$exp == (i - 1))
    if (length(idx) > 0)
      model$DLM$est[idx] <- model$DLM$est[idx] * model$Yscale / model$X[[i]]$Xscale

    for (j in i:length(model$X)) {

      if ((model$interaction > 1) || (j > i))
        model$mixNames <- c(model$mixNames,
                            paste0(model$expNames[i], "-", model$expNames[j]))

      idx <- which(model$MIX$exp1 == (i - 1) & model$MIX$exp2 == (j - 1))
      if (length(idx) > 0) {
        model$MIX$est[idx] <- model$MIX$est[idx] * model$Yscale / (model$X[[i]]$Xscale * model$X[[j]]$Xscale)
      }

    }
  }

  colnames(model$expProb) <- colnames(model$expCount) <-
    colnames(model$expInf) <- colnames(model$muExp) <- model$expNames
  if (model$interaction > 0) {
    colnames(model$mixInf) <- colnames(model$muMix) <-
      colnames(model$mixCount) <- model$mixNames
  }

  if (model$diagnostics) {
    model$treeAccept <- as.data.frame(model$treeAccept)
    colnames(model$treeAccept) <- c("tree", "step", "success", "exp", "term", "treeMhr", "mhr")
    model$treeAccept$step <- factor(model$treeAccept$step,
                                    levels = 0:3,
                                    labels = c("Grow", "Prune", "Change",
                                               "SwitchExposure"))
    model$treeAccept$exp <- factor(model$treeAccept$exp + 1,
                                   levels = 1:model$nExp,
                                   labels = model$expNames)
  }



  # Remove model and exposure data
  for (i in 1:length(model$X))
    model$X[[i]]$X <- model$X[[i]]$Tcalc <- NULL
  model$Z <- NULL

  # Change env to list
  model.out <- lapply(names(model), function(i) model[[i]])
  names(model.out) <- names(model)
  rm(model)
  gc()

  class(model.out) <- "tdlmm"
  return(model.out)
}
