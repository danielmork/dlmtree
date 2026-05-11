# Diagnostic control settings for dlmtree model fitting

Diagnostic control settings for dlmtree model fitting

## Usage

``` r
dlmtree.control.diagnose(
  subset = NULL,
  lowmem = FALSE,
  verbose = TRUE,
  save.data = TRUE,
  diagnostics = FALSE,
  initial.params = NULL
)
```

## Arguments

- subset:

  integer vector to analyze only a subset of data and exposures.

- lowmem:

  TRUE or FALSE (default): turn on memory saver for DLNM, slower
  computation time.

- verbose:

  TRUE (default) or FALSE: print output

- save.data:

  TRUE (default) or FALSE: save data used for model fitting. This must
  be set to TRUE to use shiny() function on hdlm or hdlmm

- diagnostics:

  TRUE or FALSE (default) keep model diagnostic such as the number of
  terminal nodes and acceptance ratio.

- initial.params:

  initial parameters for fixed effects model, FALSE = none (default),
  "glm" = generate using GLM, or user defined, length must equal number
  of parameters in fixed effects model.

## Value

list of control parameters for diagnostics.
