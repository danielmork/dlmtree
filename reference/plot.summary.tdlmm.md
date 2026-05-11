# Plots DLMMs for model summary of class 'tdlmm'

Method for plotting DLMMs for model summary of class 'tdlmm'. Includes
plots for marginal exposure effects as well as interactions between two
exposures.

## Usage

``` r
# S3 method for class 'summary.tdlmm'
plot(
  x,
  type = "marginal",
  exposure1 = NULL,
  exposure2 = NULL,
  time1 = c(),
  time2 = c(),
  show.cw = TRUE,
  cw.plots.only = TRUE,
  trueDLM = NULL,
  scale = NULL,
  ...
)
```

## Arguments

- x:

  an object of type 'summary.tdlmm' from summary.tdlmm() output

- type:

  plot type, 'marginal' (default)

- exposure1:

  exposure for plotting DLM

- exposure2:

  exposure paired with 'exposure1' for plotting interaction

- time1:

  plot a cross section from an interaction plot at specific time for
  'exposure1'

- time2:

  plot a cross section from an interaction plot at specific time for
  'exposure2'

- show.cw:

  indicate location of critical windows in interaction plot with red
  points

- cw.plots.only:

  show only plots with critical windows

- trueDLM:

  A vector of true effects that can be obtained from the simulated data.
  Only applicable for simulation studies

- scale:

  default = NULL, if scale is not NULL, the effects are exponentiated

- ...:

  additional plotting parameters for title and labels

## Value

A plot of distributed lag effect or interaction surface estimated with
tdlmm

## Details

plot.summary.tdlmm
