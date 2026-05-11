# diagnose

diagnose generic function for S3method

## Usage

``` r
# S3 method for class 'summary.hdlm'
diagnose(x, ...)

# S3 method for class 'summary.hdlmm'
diagnose(x, ...)

# S3 method for class 'summary.monotone'
diagnose(x, ...)

# S3 method for class 'summary.tdlm'
diagnose(x, ...)

# S3 method for class 'summary.tdlmm'
diagnose(x, ...)

# S3 method for class 'summary.tdlnm'
diagnose(x, ...)

diagnose(x, ...)
```

## Arguments

- x:

  a summary object resulting from summary() applied to an object of
  class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone'

- ...:

  not used.

## Value

shiny interface for assessing model convergence. The interface includes
tabs for MCMC diagnostics such as trace plots, density plots, and
convergence measures for distributed lag effects, DLM tree sizes, and
hyperparameters.
