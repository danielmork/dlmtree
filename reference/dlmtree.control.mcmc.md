# MCMC control settings for dlmtree model fitting

MCMC control settings for dlmtree model fitting

## Usage

``` r
dlmtree.control.mcmc(n.trees = 20, n.burn = 1000, n.iter = 2000, n.thin = 10)
```

## Arguments

- n.trees:

  integer for number of trees in ensemble.

- n.burn:

  integer for length of MCMC burn-in.

- n.iter:

  integer for number of MCMC iterations to run model after burn-in.

- n.thin:

  integer MCMC thinning factor, i.e. keep every tenth iteration.

## Value

list of MCMC control parameters.
