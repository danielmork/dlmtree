% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlmtree.R
\name{dlmtree}
\alias{dlmtree}
\title{Fit tree structured distributed lag models}
\usage{
dlmtree(
  formula,
  data,
  exposure.data,
  dlm.type = "linear",
  family = "gaussian",
  mixture = FALSE,
  het = FALSE,
  n.trees = 20,
  n.burn = 1000,
  n.iter = 2000,
  n.thin = 2,
  shrinkage = "all",
  dlmtree.params = c(0.95, 2),
  dlmtree.step.prob = c(0.25, 0.25),
  binomial.size = 1,
  formula.zi = NULL,
  tdlnm.exposure.splits = 20,
  tdlnm.time.split.prob = NULL,
  tdlnm.exposure.se = NULL,
  hdlm.modifiers = "all",
  hdlm.modifier.splits = 20,
  hdlm.modtree.params = c(0.95, 2),
  hdlm.modtree.step.prob = c(0.25, 0.25, 0.25),
  hdlm.dlmtree.type = "shared",
  hdlm.selection.prior = 0.5,
  mixture.interactions = "noself",
  mixture.prior = 1,
  monotone.gamma0 = NULL,
  monotone.sigma = NULL,
  monotone.tree.time.params = c(0.95, 2),
  monotone.tree.exp.params = c(0.95, 2),
  monotone.time.kappa = NULL,
  subset = NULL,
  lowmem = FALSE,
  verbose = TRUE,
  save.data = TRUE,
  diagnostics = FALSE,
  initial.params = NULL
)
}
\arguments{
\item{formula}{object of class formula, a symbolic description of the fixed
effect model to be fitted, e.g. y ~ a + b.}

\item{data}{data frame containing variables used in the formula.}

\item{exposure.data}{numerical matrix of exposure data with same length as data, for a mixture setting (tdlmm, hdlmm):
named list containing equally sized numerical matrices of exposure data having same length as data.}

\item{dlm.type}{dlm model specification: "linear" (default), "nonlinear", "monotone".}

\item{family}{'gaussian' for continuous response, 'logit' for binomial, 'zinb' for zero-inflated negative binomial.}

\item{mixture}{flag for mixture, set to TRUE for tdlmm and hdlmm. (default: FALSE)}

\item{het}{flag for heterogeneity, set to TRUE for hdlm and hdlmm. (default: FALSE)}

\item{n.trees}{integer for number of trees in ensemble.}

\item{n.burn}{integer for length of MCMC burn-in.}

\item{n.iter}{integer for number of MCMC iterations to run model after burn-in.}

\item{n.thin}{integer MCMC thinning factor, i.e. keep every tenth iteration.}

\item{shrinkage}{character "all" (default), "trees", "exposures", "none",
turns on horseshoe-like shrinkage priors for different parts of model.}

\item{dlmtree.params}{numerical vector of alpha and beta hyperparameters
controlling dlm tree depth. (default: alpha = 0.95, beta = 2)}

\item{dlmtree.step.prob}{numerical vector for probability of each step for dlm tree updates: 1) grow/prune,
2) change, 3) switch exposure. (default: c(0.25, 0.25, 0.25))}

\item{binomial.size}{integer type scalar (if all equal, default: 1) or
vector defining binomial size for 'logit' family.}

\item{formula.zi}{(only applies to family = 'zinb') object of class formula, a symbolic description of the fixed effect of
zero-inflated (ZI) model to be fitted, e.g. y ~ a + b. This only applies to ZINB where covariates for
ZI model are different from NB model. This is set to the argument 'formula' by default.}

\item{tdlnm.exposure.splits}{scalar indicating the number of splits (divided
evenly across quantiles of the exposure data) or list with two components:
'type' = 'values' or 'quantiles', and 'split.vals' = a numerical
vector indicating the corresponding exposure values or quantiles for splits.}

\item{tdlnm.time.split.prob}{probability vector of a spliting probabilities for time lags. (default: uniform probabilities)}

\item{tdlnm.exposure.se}{numerical matrix of exposure standard errors with same
size as exposure.data or a scalar smoothing factor representing a uniform
smoothing factor applied to each exposure measurement. (default: sd(exposure.data)/2)}

\item{hdlm.modifiers}{string vector containing desired modifiers to be included in a modifier tree.
The strings in the vector must match the names of the columns of the data.
By default, a modifier tree considers all covariates in the formula as modifiers unless stated otherwise.}

\item{hdlm.modifier.splits}{integer value to determine the possible number of splitting points that will be used for a modifier tree.}

\item{hdlm.modtree.params}{numerical vector of alpha and beta hyperparameters
controlling modifier tree depth. (default: alpha = 0.95, beta = 2)}

\item{hdlm.modtree.step.prob}{numerical vector for probability of each step for modifier tree updates: 1) grow, 2) prune,
3) change. (default: c(0.25, 0.25, 0.25))}

\item{hdlm.dlmtree.type}{specification of dlmtree type for HDLM: shared (default) or nested.}

\item{hdlm.selection.prior}{scalar hyperparameter for sparsity of modifiers. Must be between 0.5 and 1.
Smaller value corresponds to increased sparsity of modifiers.}

\item{mixture.interactions}{'noself' (default) which estimates interactions only between two
different exposures, 'all' which also allows interactions within the same exposure, or 'none'
which eliminates all interactions and estimates only main effects of each exposure.}

\item{mixture.prior}{positive scalar hyperparameter for sparsity of exposures. (default: 1)}

\item{monotone.gamma0}{vector (with length equal to number of lags) of means for logit-transformed prior probability of split at each lag;
e.g., gamma_0l = 0 implies mean prior probability of split at lag l = 0.5.}

\item{monotone.sigma}{symmetric matrix (usually with only diagonal elements) corresponding to gamma_0 to define variances on prior probability of split;
e.g., gamma_0l = 0 with lth diagonal element of sigma=2.701 implies that 95\% of the time the prior probability of split is between 0.005 and 0.995,
as a second example setting gamma_0l=4.119 and the corresponding diagonal element of sigma=0.599 implies that 95\% of the time the prior probability of a split is between 0.8 and 0.99.}

\item{monotone.tree.time.params}{numerical vector of hyperparameters for monotone time tree.}

\item{monotone.tree.exp.params}{numerical vector of hyperparameters for monotone exposure tree.}

\item{monotone.time.kappa}{scaling factor in dirichlet prior that goes alongside \code{tdlnm.time.split.prob} to
control the amount of prior information given to the model for deciding probabilities of splits between adjacent lags.}

\item{subset}{integer vector to analyze only a subset of data and exposures.}

\item{lowmem}{TRUE or FALSE (default): turn on memory saver for DLNM, slower computation time.}

\item{verbose}{TRUE (default) or FALSE: print output}

\item{save.data}{TRUE (default) or FALSE: save data used for model fitting. This must be set to TRUE to use shiny() function on hdlm or hdlmm}

\item{diagnostics}{TRUE or FALSE (default) keep model diagnostic such as the number of
terminal nodes and acceptance ratio.}

\item{initial.params}{initial parameters for fixed effects model, FALSE = none (default),
"glm" = generate using GLM, or user defined, length must equal number of parameters in fixed effects model.}
}
\value{
Object of one of the classes: tdlm, tdlmm, tdlnm, hdlm, hdlmm
}
\description{
The 'dlmtree' function accommodates various response variable types, including continuous, binary, and zero-inflated count values.
The function is designed to handle both single exposure and exposure mixtures. For a single exposure, users are offered options to model non-linear effects (tdlnm),
linear effects (tdlm), or heterogeneous subgroup/individualized effects (hdlm). In the case of exposure mixtures, the function supports
lagged interactions (tdlmm), and heterogeneous subgroup/individualized effects (hdlmm) allowing for a comprehensive exploration of mixture exposure heterogeneity.
Additionally, users can fine-tune parameters to impose effect shrinkage and perform exposure selection, enhancing the adaptability and precision of the modeling process.
For more detailed documentation, visit: \href{https://danielmork.github.io/dlmtree/}{dlmtree website}.
}
\details{
dlmtree

Model is recommended to be run for at minimum 5000 burn-in
iterations followed by 15000 sampling iterations with a thinning factor of 5.
Convergence can be checked by re-running the model and validating consistency
of results. Examples are provided below for the syntax for running different types of models.
For more examples, visit: \href{https://danielmork.github.io/dlmtree/}{dlmtree website}.
}
\examples{
\donttest{

  # The first three examples are for one lagged exposure


  # treed distributed lag model (TDLM)
  # binary outcome with logit link

  D <- sim.tdlmm(sim = "A", mean.p = 0.5, n = 1000)
  tdlm.fit <- dlmtree(y ~ .,
                      data = D$dat,
                      exposure.data = D$exposures[[1]],
                      dlm.type = "linear",
                      family = "logit",
                      binomial.size = 1)

  # summarize results
  tdlm.sum <- summary(tdlm.fit)
  tdlm.sum

  # plot results
  plot(tdlm.sum)



  # Treed distributed lag nonlinear model (TDLNM)
  # Gaussian regression model
  D <- sim.tdlnm(sim = "A", error.to.signal = 1)
  tdlnm.fit <- dlmtree(formula = y ~ .,
                       data = D$dat,
                       exposure.data = D$exposures,
                       dlm.type = "nonlinear",
                       family = "gaussian")

  # summarize results
  tdlnm.sum <- summary(tdlnm.fit)
  tdlnm.sum

  # plot results
  plot(tdlnm.sum)



  # Heterogenious TDLM (HDLM), similar to first example but with heterogenious exposure response
  D <- sim.hdlmm(sim = "B", n = 1000)
  hdlm.fit <- dlmtree(y ~ .,
                      data = D$dat,
                      exposure.data = D$exposures,
                      dlm.type = "linear",
                      family = "gaussian",
                      het = TRUE)

  # summarize results
  hdlm.sum <- summary(hdlm.fit)
  hdlm.sum

  # shiny app for HDLM
  if (interactive()) {
    shiny(hdlm.fit)
  }



  # The next two examples are for a mixture (or multivariate) exposure


  # Treed distributed lag mixture model (TDLMM)
  # Model for mixutre (or multivariate) lagged exposures
  # with a homogenious exposure-time-response function
  D <- sim.tdlmm(sim = "B", error = 25, n = 1000)
  tdlmm.fit <- dlmtree(y ~ .,
                       data = D$dat, exposure.data = D$exposures,
                       mixture.interactions = "noself",
                       dlm.type = "linear", family = "gaussian",
                       mixture = TRUE)

  # summarize results
  tdlmm.sum <- summary(tdlmm.fit)

  # plot the marginal exposure-response for one exposure
  plot(tdlmm.sum, exposure1 = "e1")

  # plot exposure-response surface
  plot(tdlmm.sum, exposure1 = "e1", exposure2 = "e2")



  # heterogenious version of TDLMM
  D <- sim.hdlmm(sim = "D", n = 1000)
  hdlmm.fit <- dlmtree(y ~ .,
                       data = D$dat,
                       exposure.data = D$exposures,
                       dlm.type = "linear",
                       family = "gaussian",
                       mixture = TRUE,
                       het = TRUE)

  # summarize results
  hdlmm.sum <- summary(hdlmm.fit)
  hdlmm.sum

  # summarize results
  if (interactive()) {
    shiny(hdlmm.fit)
  }


}
}
