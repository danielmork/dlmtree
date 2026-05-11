# Returns variety of plots for model summary of class 'monotone'

Method for returning variety of plots for model summary of class
'monotone'

## Usage

``` r
# S3 method for class 'summary.monotone'
plot(x, plot.type = "mean", val = c(), time = c(), ...)
```

## Arguments

- x:

  object of class 'summary.monotone', output of summary of 'monotone'

- plot.type:

  string indicating plot type, options are 'mean' (default) which shows
  mean exposure-time response surface, 'se', 'ci-min', 'ci-max', 'slice'
  which takes a slice of the plot at a given 'val' or 'time', 'animate'
  which creates a animation of slices of the surface plot across
  exposure values (requires package gganimate)

- val:

  exposure value for slice plot

- time:

  time value for slice plot

- ...:

  additional parameters to alter plots: 'main', 'xlab', 'ylab', 'flab'
  which sets the effect label for surface plots, 'start.time' which sets
  the first time value

## Value

A plot of distributed lag effect estimated with monotone-TDLNM

## Details

plot.summary.monotone
