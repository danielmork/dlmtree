# Returns variety of plots for model summary of class 'tdlnm'

Method for returning variety of plots for model summary of class 'tdlnm'

## Usage

``` r
# S3 method for class 'summary.tdlnm'
plot(x, plot.type = "mean", val = c(), time = c(), ...)
```

## Arguments

- x:

  object of class 'summary.tdlnm', output of summary of 'tdlnm'

- plot.type:

  string indicating plot type, options are 'mean' (default) which shows
  mean exposure-time response surface, 'cumulative' which shows the
  cumulative effects per exposure-concentration level, 'effect' which
  returns a grid of exposure concentration and lag to determine if the
  credible interval contains zero, with the direction of the effect
  indicated, 'se', 'ci-min', 'ci-max', 'slice' which takes a slice of
  the plot at a given 'val' or 'time', 'animate' which creates a
  animation of slices of the surface plot across exposure values
  (requires package gganimate)

- val:

  exposure value for slice plot

- time:

  time value for slice plot

- ...:

  additional plotting parameters for title and labels 'flab' which sets
  the effect label for surface plots, 'start.time' which sets the first
  time value

## Value

A plot of distributed lag effect estimated with tdlnm

## Details

plot.summary.tdlnm
