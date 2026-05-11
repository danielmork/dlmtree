# Hyperparameter control settings for dlmtree model fitting

Hyperparameter control settings for dlmtree model fitting

## Usage

``` r
dlmtree.control.hyper(
  shrinkage = "all",
  params = c(0.95, 2),
  step.prob = c(0.25, 0.25)
)
```

## Arguments

- shrinkage:

  character "all" (default), "trees", "exposures", "none", turns on
  horseshoe-like shrinkage priors for different parts of model.

- params:

  numerical vector of alpha and beta hyperparameters controlling dlm
  tree depth. (default: alpha = 0.95, beta = 2)

- step.prob:

  numerical vector for probability of each step for dlm tree updates: 1)
  grow/prune, 2) change, 3) switch exposure. (default: c(0.25, 0.25,
  0.25))

## Value

list of hyperparameter control parameters.
