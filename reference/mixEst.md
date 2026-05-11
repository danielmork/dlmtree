# Calculates the lagged interaction effects with MIX matrix for linear models.

Calculates the lagged interaction effects with MIX matrix for linear
models.

## Usage

``` r
mixEst(dlm, nlags, nsamp)
```

## Arguments

- dlm:

  A numeric matrix containing the model fit information

- nlags:

  total number of lags

- nsamp:

  number of mcmc iterations

## Value

A cube object of interaction effect x lag x mcmc
