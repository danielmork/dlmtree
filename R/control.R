#' MCMC control settings for dlmtree model fitting
#' 
#' @param n.trees integer for number of trees in ensemble.
#' @param n.burn integer for length of MCMC burn-in.
#' @param n.iter integer for number of MCMC iterations to run model after burn-in.
#' @param n.thin integer MCMC thinning factor, i.e. keep every tenth iteration.
#' 
#' @returns A list of MCMC control parameters.
#' @export
ctrl.mcmc <- function(
    n.trees = 20,
    n.burn = 1000,
    n.iter = 2000,
    n.thin = 10
) {
 list(
    n.trees = n.trees,
    n.burn = n.burn,
    n.iter = n.iter,
    n.thin = n.thin
  )
}

#' Hyperparameter control settings for dlmtree model fitting
#' 
#' @param shrinkage character "all" (default), "trees", "exposures", "none",
#' turns on horseshoe-like shrinkage priors for different parts of model.
#' @param dlmtree.params numerical vector of alpha and beta hyperparameters
#' controlling dlm tree depth. (default: alpha = 0.95, beta = 2)
#' @param dlmtree.step.prob numerical vector for probability of each step for dlm tree updates: 1) grow/prune,
#' 2) change, 3) switch exposure. (default: c(0.25, 0.25, 0.25))
#' 
#' @returns A list of hyperparameter control parameters.
#' @export
ctrl.hyperparam <- function(
    shrinkage = "all",
    dlmtree.params = c(.95, 2),
    dlmtree.step.prob = c(.25, .25)
){
  list(
    shrinkage = shrinkage,
    dlmtree.params = dlmtree.params,
    dlmtree.step.prob = dlmtree.step.prob
  )
}

#' Family control settings for dlmtree model fitting
#' 
#' @param binomial.size integer type scalar (if all equal, default: 1) or vector defining binomial size for 'logit' family.
#' @param formula.zi (only applies to family = 'zinb') object of class formula, a symbolic description of the fixed effect of
#' zero-inflated (ZI) model to be fitted, e.g. y ~ a + b. This only applies to ZINB where covariates for
#' ZI model are different from NB model. This is set to the argument 'formula' by default.
#' 
#' @returns A list of family control parameters.
#' @export
ctrl.family <- function(
    binomial.size = 1,
    formula.zi = NULL
) {
  list(
    binomial.size = binomial.size, 
    formula.zi = formula.zi)
}

#' Control settings for dlmtree model fitting, when used for TDLNM
#' 
#' @param tdlnm.exposure.splits scalar indicating the number of splits (divided
#' evenly across quantiles of the exposure data) or list with two components:
#' 'type' = 'values' or 'quantiles', and 'split.vals' = a numerical
#' vector indicating the corresponding exposure values or quantiles for splits.
#' @param tdlnm.exposure.se numerical matrix of exposure standard errors with same
#' size as exposure.data or a scalar smoothing factor representing a uniform
#' smoothing factor applied to each exposure measurement. (default: sd(exposure.data)/2)
#' @param tdlnm.time.split.prob probability vector of a spliting probabilities for time lags. (default: uniform probabilities)
#' 
#' @returns A list of TDLNM control parameters.
#' @export
ctrl.tdlnm <- function(
    tdlnm.exposure.splits = 20,
    tdlnm.time.split.prob = NULL,
    tdlnm.exposure.se = NULL
) {
  list(
    tdlnm.exposure.splits = tdlnm.exposure.splits,
    tdlnm.time.split.prob = tdlnm.time.split.prob,
    tdlnm.exposure.se = tdlnm.exposure.se
  )
}


#' Control settings for dlmtree model fitting, when used for heterogeneous models
#' 
#' @param hdlm.modifiers string vector containing desired modifiers to be included in a modifier tree.
#' The strings in the vector must match the names of the columns of the data. 
#' By default, a modifier tree considers all covariates in the formula as modifiers unless stated otherwise.
#' @param hdlm.modifier.splits integer value to determine the possible number of splitting points that will be used for a modifier tree.
#' @param hdlm.modtree.params numerical vector of alpha and beta hyperparameters
#' controlling modifier tree depth. (default: alpha = 0.95, beta = 2)
#' @param hdlm.modtree.step.prob numerical vector for probability of each step for modifier tree updates: 1) grow, 2) prune,
#' 3) change. (default: c(0.25, 0.25, 0.25))
#' @param hdlm.dlmtree.type specification of dlmtree type for HDLM: shared (default) or nested.
#' @param hdlm.selection.prior scalar hyperparameter for sparsity of modifiers. Must be between 0.5 and 1. 
#' Smaller value corresponds to increased sparsity of modifiers.
#' 
#' @returns A list of control parameters for heterogeneous models.
#' @export
ctrl.het <- function(
    hdlm.modifiers = "all",
    hdlm.modifier.splits = 20,
    hdlm.modtree.params = c(.95, 2),
    hdlm.modtree.step.prob = c(.25, .25, .25),
    hdlm.dlmtree.type = "shared",
    hdlm.selection.prior = 0.5
) {
  list(
    hdlm.modifiers = hdlm.modifiers,
    hdlm.modifier.splits = hdlm.modifier.splits,
    hdlm.modtree.params = hdlm.modtree.params,
    hdlm.modtree.step.prob = hdlm.modtree.step.prob,
    hdlm.dlmtree.type = hdlm.dlmtree.type,
    hdlm.selection.prior = hdlm.selection.prior
  )
}

#' Control settings for dlmtree model fitting, when used for mixture models
#' 
#' @param mixture.interactions 'noself' (default) which estimates interactions only between two 
#' different exposures, 'all' which also allows interactions within the same exposure, or 'none' 
#' which eliminates all interactions and estimates only main effects of each exposure.
#' @param mixture.prior positive scalar hyperparameter for sparsity of exposures. (default: 1)
#' 
#' @returns A list of mixture control parameters.
#' @export
ctrl.mixture <- function(
    mixture.interactions = "noself",
    mixture.prior = 1
) {
  list(
    mixture.interactions = mixture.interactions,
    mixture.prior = mixture.prior)
}

#' Control settings for dlmtree model fitting, when used for monotone model
#' 
#' @param monotone.gamma0 vector (with length equal to number of lags) of means for logit-transformed prior probability of split at each lag; e.g., gamma_0l = 0 implies mean prior probability of split at lag l = 0.5.
#' @param monotone.sigma symmetric matrix (usually with only diagonal elements) corresponding to gamma_0 to define variances on prior probability of split; 
#' e.g., gamma_0l = 0 with lth diagonal element of sigma=2.701 implies that 95\% of the time the prior probability of split is between 0.005 and 0.995, 
#' as a second example setting gamma_0l=4.119 and the corresponding diagonal element of sigma=0.599 implies that 95\% of the time the prior probability of a split is between 0.8 and 0.99.
#' @param monotone.tree.time.params numerical vector of hyperparameters for monotone time tree.
#' @param monotone.tree.exp.params numerical vector of hyperparameters for monotone exposure tree.
#' @param monotone.time.kappa scaling factor in dirichlet prior that goes alongside `tdlnm.time.split.prob` to control the amount of prior information given to the model for deciding probabilities of splits between adjacent lags.
#' 
#' @returns A list of control parameters for monotone model.
#' @export
ctrl.monotone <- function(
    monotone.gamma0 = NULL,
    monotone.sigma = NULL,
    monotone.tree.time.params = c(.95, 2),
    monotone.tree.exp.params = c(.95, 2),
    monotone.time.kappa = NULL
) {
  list(
    monotone.gamma0 = monotone.gamma0,
    monotone.sigma = monotone.sigma,
    monotone.tree.time.params = monotone.tree.time.params,
    monotone.tree.exp.params = monotone.tree.exp.params,
    monotone.time.kappa = monotone.time.kappa
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
#' @returns A list of control parameters for diagnostics.
#' @export
ctrl.diagnose <- function(
    subset = NULL,
    lowmem = FALSE,
    verbose = TRUE,
    save.data = TRUE,
    diagnostics = FALSE,
    initial.params = NULL
) {
  list(
    subset = subset,
    lowmem = lowmem,
    verbose = verbose,
    save.data = save.data,
    diagnostics = diagnostics,
    initial.params = initial.params
  )
}
