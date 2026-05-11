# summary

summary generic function for S3method

## Usage

``` r
summary(x, conf.level = 0.95, ...)

# S3 method for class 'hdlm'
summary(x, conf.level = 0.95, mcmc = FALSE, ...)

# S3 method for class 'hdlmm'
summary(x, conf.level = 0.95, mcmc = FALSE, ...)

# S3 method for class 'monotone'
summary(
  x,
  conf.level = 0.95,
  pred.at = NULL,
  cenval = 0,
  exposure.se = NULL,
  mcmc = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'tdlm'
summary(x, conf.level = 0.95, mcmc = FALSE, ...)

# S3 method for class 'tdlmm'
summary(
  x,
  conf.level = 0.95,
  marginalize = "mean",
  log10BF.crit = 0.5,
  mcmc = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'tdlnm'
summary(
  x,
  conf.level = 0.95,
  pred.at = NULL,
  cenval = 0,
  exposure.se = NULL,
  mcmc = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  an object of class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm',
  'monotone'

- conf.level:

  confidence level for computation of credible intervals

- ...:

  additional parameters

- mcmc:

  keep all mcmc iterations (large memory requirement)

- pred.at:

  numerical vector of exposure values to make predictions for at each
  time period

- cenval:

  scalar exposure value that acts as a reference point for predictions
  at all other exposure values

- exposure.se:

  scalar smoothing factor, if different from model

- verbose:

  show progress in console

- marginalize:

  value(s) for calculating marginal DLMs, defaults to "mean", can also
  specify a percentile from 1-99 for all other exposures, or a named
  vector with specific values for each exposure

- log10BF.crit:

  Bayes Factor criteria for selecting exposures and interactions, such
  that log10(BayesFactor) \> x. Default = 0.5.

## Value

list of summary outputs of the model fit
