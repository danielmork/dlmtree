# Determines split points for continuous modifiers

Method for determining split points for continuous modifiers

## Usage

``` r
splitpoints(object, var, round = NULL)
```

## Arguments

- object:

  An object of class 'hdlm', 'hdlmm'

- var:

  The name of a continuous variable for which the split points will be
  reported

- round:

  The number of decimal places to round the variable (var) to. No
  rounding occurs if round=NULL (default) For positive integer values of
  round, the variable will be rounded and split points will be reported
  at the resulting level

## Value

A data frame with split points and the probability that a split point
was \>= that split point value

## Details

splitpoints

## Examples

``` r
# \donttest{
# Split points with HDLM 
D <- sim.hdlmm(sim = "B", n = 1000)
fit <- dlmtree(y ~ ., 
               data = D$dat,
               exposure.data = D$exposures,
               dlm.type = "linear",
               family = "gaussian",
               het = TRUE)
#> Preparing data...
#> 
#> Running shared HDLM:
#> Burn-in % complete 
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> MCMC iterations (est time: 8 seconds)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
splitpoints(fit, var = "mod_num", round = 2)
#>    location  proportion
#> 1     -1.75 0.041121495
#> 2     -1.33 0.034579439
#> 3     -1.06 0.024766355
#> 4     -0.89 0.025233645
#> 5     -0.75 0.071495327
#> 6     -0.63 0.107009346
#> 7     -0.48 0.050000000
#> 8     -0.35 0.035514019
#> 9     -0.22 0.040654206
#> 10    -0.10 0.167289720
#> 11     0.02 0.128504673
#> 12     0.17 0.079439252
#> 13     0.30 0.035514019
#> 14     0.42 0.017289720
#> 15     0.51 0.021028037
#> 16     0.65 0.022429907
#> 17     0.78 0.018691589
#> 18     0.97 0.036915888
#> 19     1.19 0.039719626
#> 20     1.57 0.002803738
splitpoints(fit, var = "mod_scale", round = 2)
#>    location proportion
#> 1      0.05 0.03118416
#> 2      0.09 0.03076275
#> 3      0.16 0.09734513
#> 4      0.20 0.03329119
#> 5      0.25 0.02907712
#> 6      0.30 0.01980615
#> 7      0.36 0.03581964
#> 8      0.41 0.10408765
#> 9      0.46 0.02907712
#> 10     0.50 0.02275601
#> 11     0.55 0.01980615
#> 12     0.59 0.06110409
#> 13     0.63 0.02823430
#> 14     0.68 0.03329119
#> 15     0.72 0.03286979
#> 16     0.78 0.05731142
#> 17     0.82 0.05520438
#> 18     0.86 0.06068268
#> 19     0.91 0.16982722
#> 20     0.97 0.04846186
# }
```
