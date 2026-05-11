# Calculates the distributed lag effect with DLM matrix for non-linear models.

Calculates the distributed lag effect with DLM matrix for non-linear
models.

## Usage

``` r
dlnmPLEst(dlnm, predAt, nlags, nsamp, center)
```

## Arguments

- dlnm:

  A numeric matrix containing the model fit information

- predAt:

  Number of splits in the model

- nlags:

  total number of lags

- nsamp:

  number of mcmc iterations

- center:

  center parameter

## Value

A cube object of lag effect x lag x mcmc
