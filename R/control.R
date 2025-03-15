#' MCMC control settings for dlmtree model fitting
#' 
#' @param n.trees integer for number of trees in ensemble.
#' @param n.burn integer for length of MCMC burn-in.
#' @param n.iter integer for number of MCMC iterations to run model after burn-in.
#' @param n.thin integer MCMC thinning factor, i.e. keep every tenth iteration.
#' 
#' @returns list of MCMC control parameters.
#' @export
dlmtree.control.mcmc <- function(
    n.trees = 20,
    n.burn  = 1000,
    n.iter  = 2000,
    n.thin  = 10
) {
 list(
    n.trees = n.trees,
    n.burn  = n.burn,
    n.iter  = n.iter,
    n.thin  = n.thin
  )
}

#' Hyperparameter control settings for dlmtree model fitting
#' 
#' @param shrinkage character "all" (default), "trees", "exposures", "none",
#' turns on horseshoe-like shrinkage priors for different parts of model.
#' @param params numerical vector of alpha and beta hyperparameters
#' controlling dlm tree depth. (default: alpha = 0.95, beta = 2)
#' @param step.prob numerical vector for probability of each step for dlm tree updates: 1) grow/prune,
#' 2) change, 3) switch exposure. (default: c(0.25, 0.25, 0.25))
#' 
#' @returns list of hyperparameter control parameters.
#' @export
dlmtree.control.hyper <- function(
    shrinkage = "all",
    params    = c(.95, 2),
    step.prob = c(.25, .25)
){
  list(
    shrinkage = shrinkage,
    params    = params,
    step.prob = step.prob
  )
}

#' Family control settings for dlmtree model fitting
#' 
#' @param binomial.size integer type scalar (if all equal, default: 1) or vector defining binomial size for 'logit' family.
#' @param formula.zi (only applies to family = 'zinb') object of class formula, a symbolic description of the fixed effect of
#' zero-inflated (ZI) model to be fitted, e.g. y ~ a + b. This only applies to ZINB where covariates for
#' ZI model are different from NB model. This is set to the argument 'formula' by default.
#' 
#' @returns list of family control parameters.
#' @export
dlmtree.control.family <- function(
    binomial.size = 1,
    formula.zi    = NULL
) {
  list(
    binomial.size = binomial.size, 
    formula.zi    = formula.zi)
}

#' Control settings for dlmtree model fitting, when used for TDLNM
#' 
#' @param exposure.splits scalar indicating the number of splits (divided
#' evenly across quantiles of the exposure data) or list with two components:
#' 'type' = 'values' or 'quantiles', and 'split.vals' = a numerical
#' vector indicating the corresponding exposure values or quantiles for splits.
#' @param exposure.se numerical matrix of exposure standard errors with same
#' size as exposure.data or a scalar smoothing factor representing a uniform
#' smoothing factor applied to each exposure measurement. (default: sd(exposure.data)/2)
#' @param time.split.prob probability vector of a spliting probabilities for time lags. (default: uniform probabilities)
#' 
#' @returns list of TDLNM control parameters.
#' @export
dlmtree.control.tdlnm <- function(
    exposure.splits = 20,
    time.split.prob = NULL,
    exposure.se     = NULL
) {
  list(
    exposure.splits = exposure.splits,
    time.split.prob = time.split.prob,
    exposure.se     = exposure.se
  )
}


#' Control settings for dlmtree model fitting, when used for heterogeneous models
#' 
#' @param modifiers string vector containing desired modifiers to be included in a modifier tree.
#' The strings in the vector must match the names of the columns of the data. 
#' By default, a modifier tree considers all covariates in the formula as modifiers unless stated otherwise.
#' @param modifier.splits integer value to determine the possible number of splitting points that will be used for a modifier tree.
#' @param modtree.params numerical vector of alpha and beta hyperparameters
#' controlling modifier tree depth. (default: alpha = 0.95, beta = 2)
#' @param modtree.step.prob numerical vector for probability of each step for modifier tree updates: 1) grow, 2) prune,
#' 3) change. (default: c(0.25, 0.25, 0.25))
#' @param dlmtree.type specification of dlmtree type for HDLM: shared (default) or nested.
#' @param selection.prior scalar hyperparameter for sparsity of modifiers. Must be between 0.5 and 1. 
#' Smaller value corresponds to increased sparsity of modifiers.
#' 
#' @returns list of control parameters for heterogeneous models.
#' @export
dlmtree.control.het <- function(
    modifiers         = "all",
    modifier.splits   = 20,
    modtree.params    = c(.95, 2),
    modtree.step.prob = c(.25, .25, .25),
    dlmtree.type      = "shared",
    selection.prior   = 0.5
) {
  list(
    modifiers         = modifiers,
    modifier.splits   = modifier.splits,
    modtree.params    = modtree.params,
    modtree.step.prob = modtree.step.prob,
    dlmtree.type      = dlmtree.type,
    selection.prior   = selection.prior
  )
}

#' Control settings for dlmtree model fitting, when used for mixture models
#' 
#' @param interactions 'noself' (default) which estimates interactions only between two 
#' different exposures, 'all' which also allows interactions within the same exposure, or 'none' 
#' which eliminates all interactions and estimates only main effects of each exposure.
#' @param sparsity.prior positive scalar hyperparameter for sparsity of exposures. (default: 1)
#' 
#' @returns list of mixture control parameters.
#' @export
dlmtree.control.mix <- function(
    interactions   = "noself",
    sparsity.prior = 1
) {
  list(
    interactions   = interactions,
    sparsity.prior = sparsity.prior)
}

#' Control settings for dlmtree model fitting, when used for monotone model
#' 
#' @param gamma0 vector (with length equal to number of lags) of means for logit-transformed prior probability of split at each lag; e.g., gamma_0l = 0 implies mean prior probability of split at lag l = 0.5.
#' @param sigma symmetric matrix (usually with only diagonal elements) corresponding to gamma_0 to define variances on prior probability of split; 
#' e.g., gamma_0l = 0 with lth diagonal element of sigma=2.701 implies that 95\% of the time the prior probability of split is between 0.005 and 0.995, 
#' as a second example setting gamma_0l=4.119 and the corresponding diagonal element of sigma=0.599 implies that 95\% of the time the prior probability of a split is between 0.8 and 0.99.
#' @param tree.time.params numerical vector of hyperparameters for monotone time tree.
#' @param tree.exp.params numerical vector of hyperparameters for monotone exposure tree.
#' @param time.kappa scaling factor in dirichlet prior that goes alongside `time.split.prob` to control the amount of prior information given to the model for deciding probabilities of splits between adjacent lags.
#' 
#' @returns list of control parameters for monotone model.
#' @export
dlmtree.control.monotone <- function(
    gamma0           = NULL,
    sigma            = NULL,
    tree.time.params = c(.95, 2),
    tree.exp.params  = c(.95, 2),
    time.kappa       = NULL
) {
  list(
    gamma0           = gamma0,
    sigma            = sigma,
    tree.time.params = tree.time.params,
    tree.exp.params  = tree.exp.params,
    time.kappa       = time.kappa
  )
}

#' Diagnostic control settings for dlmtree model fitting
#' 
#' @param subset integer vector to analyze only a subset of data and exposures.
#' @param lowmem TRUE or FALSE (default): turn on memory saver for DLNM, slower computation time.
#' @param verbose TRUE (default) or FALSE: print output
#' @param save.data TRUE (default) or FALSE: save data used for model fitting. This must be set to TRUE to use shiny() function on hdlm or hdlmm
#' @param diagnostics TRUE or FALSE (default) keep model diagnostic such as the number of
#' terminal nodes and acceptance ratio.
#' @param initial.params initial parameters for fixed effects model, FALSE = none (default), 
#' "glm" = generate using GLM, or user defined, length must equal number of parameters in fixed effects model.
#' 
#' @returns list of control parameters for diagnostics.
#' @export
dlmtree.control.diagnose <- function(
    subset         = NULL,
    lowmem         = FALSE,
    verbose        = TRUE,
    save.data      = TRUE,
    diagnostics    = FALSE,
    initial.params = NULL
) {
  list(
    subset         = subset,
    lowmem         = lowmem,
    verbose        = verbose,
    save.data      = save.data,
    diagnostics    = diagnostics,
    initial.params = initial.params
  )
}
