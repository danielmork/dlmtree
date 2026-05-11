# Control settings for dlmtree model fitting, when used for mixture models

Control settings for dlmtree model fitting, when used for mixture models

## Usage

``` r
dlmtree.control.mix(interactions = "noself", sparsity.prior = 1)
```

## Arguments

- interactions:

  'noself' (default) which estimates interactions only between two
  different exposures, 'all' which also allows interactions within the
  same exposure, or 'none' which eliminates all interactions and
  estimates only main effects of each exposure.

- sparsity.prior:

  positive scalar hyperparameter for sparsity of exposures. (default: 1)

## Value

list of mixture control parameters.
