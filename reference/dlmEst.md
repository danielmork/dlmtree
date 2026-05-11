# Calculates the distributed lag effect with DLM matrix for linear models.

Calculates the distributed lag effect with DLM matrix for linear models.

## Usage

``` r
dlmEst(dlm, nlags, nsamp)
```

## Arguments

- dlm:

  A numeric matrix containing the model fit information

- nlags:

  total number of lags

- nsamp:

  number of mcmc iterations

## Value

A cube object of lag effect x lag x mcmc
