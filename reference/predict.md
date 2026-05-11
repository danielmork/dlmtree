# predict

predict generic function for S3method

## Usage

``` r
predict(
  x,
  new.data,
  new.exposure.data,
  ci.level = 0.95,
  type = "response",
  outcome = NULL,
  fixed.idx = list(),
  est.dlm = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'hdlm'
predict(
  x,
  new.data,
  new.exposure.data,
  ci.level = 0.95,
  type = "response",
  outcome = NULL,
  fixed.idx = list(),
  est.dlm = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'hdlmm'
predict(
  x,
  new.data,
  new.exposure.data,
  ci.level = 0.95,
  type = "response",
  outcome = NULL,
  fixed.idx = list(),
  est.dlm = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  fitted dlmtree model with class 'hdlm', 'hdlmm'

- new.data:

  new data frame which contains the same covariates and modifiers used
  to fit the model

- new.exposure.data:

  new data frame/list which contains the same length of exposure lags
  used to fit the model

- ci.level:

  credible interval level for posterior predictive distribution

- type:

  type of prediction: "response" (default) or "waic". "waic" must be
  specified with \`outcome\` parameter

- outcome:

  outcome required for WAIC calculation

- fixed.idx:

  fixed index

- est.dlm:

  flag for estimating dlm effect

- verbose:

  TRUE (default) or FALSE: print output

- ...:

  not used

## Value

list with the following elements:

- ztg:

  posterior predictive mean of fixed effect

- ztg.lims:

  lower/upper bound of posterior predictive distribution of fixed effect

- dlmest:

  estimated exposure effect

- dlmest.lower:

  lower bound of estimated exposure effect

- dlmest.upper:

  upper bound of estimated exposure effect

- fhat:

  posterior predictive mean of exposure effect

- fhat.lims:

  lower/upper bound of posterior predictive distribution of exposure
  effect

- y:

  posterior predictive mean

- y.lims:

  lower/upper bound of posterior predictive distribution
