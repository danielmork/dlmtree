# Calculates subgroup-specific lag effects for heterogeneous models

Method for calculating subgroup-specific lag effects for heterogeneous
models: HDLM, HDLMM

## Usage

``` r
estDLM(
  object,
  new.data,
  group.index,
  conf.level = 0.95,
  exposure = NULL,
  keep.mcmc = FALSE,
  mem.safe = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  object of a model fit. Must be 'hdlm' or 'hdlmm'

- new.data:

  data frame with new observations with the same number of modifiers

- group.index:

  list of index (row numbers) for subgroup specification

- conf.level:

  confidence level for credible interval of effects

- exposure:

  exposure of interest for 'hdlmm' method

- keep.mcmc:

  store mcmc in the output

- mem.safe:

  boolean memory parameter for rule index

- verbose:

  TRUE (default) or FALSE: print output

## Value

A list with the following components:

- conf.level:

  Specified confidence level

- mod:

  a list of modifers with a vector of values from the model

- n:

  Number of observation per specified subgroup

- groupIndex:

  list of index (row numbers) for specified subgroup

- dlmMean:

  distributed lag effects per subgroups

- dlmCI:

  credible intervals for distributed lag effects per subgroups

- dlmCum:

  cumulative effects per subgroups

- dlFunction:

  type of DLM class

- plotData:

  data frame built for easier visualization of distributed lag effects
  for each subgroup (facet)

## Details

estDLM
