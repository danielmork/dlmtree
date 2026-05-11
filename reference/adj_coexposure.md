# Adjusting for expected changes in co-exposure (TDLMM)

Estimates the marginal effects of an exposure while accounting for
expected changes in co-occurring exposures at the same time point.
Values of co-occurring exposures are modeled nonlinearly using a spline
model with predictions made at the lower an upper values for the
exposure of interest.

## Usage

``` r
adj_coexposure(
  exposure.data,
  object,
  contrast_perc = c(0.25, 0.75),
  contrast_exp = list(),
  conf.level = 0.95,
  keep.mcmc = FALSE,
  verbose = TRUE
)
```

## Arguments

- exposure.data:

  Named list of exposure matrices used as input to TDLMM.

- object:

  Model output for TDLMM from dlmtree() function.

- contrast_perc:

  2-length vector of percentiles or named list corresponding to lower
  and upper exposure percentiles of interest. Names must equal list
  names in 'exposure.data'.

- contrast_exp:

  Named list consisting lower and upper exposure values. This takes
  precedence over contrast_perc if both inputs are used.

- conf.level:

  Confidence level used for estimating credible intervals. Default is
  0.95.

- keep.mcmc:

  If TRUE, return posterior samples.

- verbose:

  TRUE (default) or FALSE: print output

## Value

A list with the following components (or posterior samples if keep.mcmc
= TRUE):

- Name:

  vector of exposure names

- Time:

  integer vector of lags

- Effect:

  posterior mean of marginal effects

- SE:

  standard error of the estimate

- Lower:

  lower bound of credible interval of the marginal effect estimate

- Upper:

  upper bound of credible interval of the marginal effect estimate

- cEffect:

  cumulative marginal effects

- cLower:

  lower bound of credible interval of the cumulative marginal effect

- cUpper:

  upper bound of credible interval of the cumulative marginal effect

- CW:

  boolean vector indicating critical window

## Details

adj_coexposure
