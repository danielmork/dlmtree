# Control settings for dlmtree model fitting, when used for monotone model

Control settings for dlmtree model fitting, when used for monotone model

## Usage

``` r
dlmtree.control.monotone(
  gamma0 = NULL,
  sigma = NULL,
  tree.time.params = c(0.95, 2),
  tree.exp.params = c(0.95, 2),
  time.kappa = NULL
)
```

## Arguments

- gamma0:

  vector (with length equal to number of lags) of means for
  logit-transformed prior probability of split at each lag; e.g.,
  gamma_0l = 0 implies mean prior probability of split at lag l = 0.5.

- sigma:

  symmetric matrix (usually with only diagonal elements) corresponding
  to gamma_0 to define variances on prior probability of split; e.g.,
  gamma_0l = 0 with lth diagonal element of sigma=2.701 implies that 95%
  of the time the prior probability of split is between 0.005 and 0.995,
  as a second example setting gamma_0l=4.119 and the corresponding
  diagonal element of sigma=0.599 implies that 95% of the time the prior
  probability of a split is between 0.8 and 0.99.

- tree.time.params:

  numerical vector of hyperparameters for monotone time tree.

- tree.exp.params:

  numerical vector of hyperparameters for monotone exposure tree.

- time.kappa:

  scaling factor in dirichlet prior that goes alongside
  \`time.split.prob\` to control the amount of prior information given
  to the model for deciding probabilities of splits between adjacent
  lags.

## Value

list of control parameters for monotone model.
