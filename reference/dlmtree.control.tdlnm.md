# Control settings for dlmtree model fitting, when used for TDLNM

Control settings for dlmtree model fitting, when used for TDLNM

## Usage

``` r
dlmtree.control.tdlnm(
  exposure.splits = 20,
  time.split.prob = NULL,
  exposure.se = NULL
)
```

## Arguments

- exposure.splits:

  scalar indicating the number of splits (divided evenly across
  quantiles of the exposure data) or list with two components: 'type' =
  'values' or 'quantiles', and 'split.vals' = a numerical vector
  indicating the corresponding exposure values or quantiles for splits.

- time.split.prob:

  probability vector of a spliting probabilities for time lags.
  (default: uniform probabilities)

- exposure.se:

  numerical matrix of exposure standard errors with same size as
  exposure.data or a scalar smoothing factor representing a uniform
  smoothing factor applied to each exposure measurement. (default:
  sd(exposure.data)/2)

## Value

list of TDLNM control parameters.
