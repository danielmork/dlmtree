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
#> MCMC iterations (est time: 6 seconds)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
splitpoints(fit, var = "mod_num", round = 2)
#>    location  proportion
#> 1     -1.69 0.043936093
#> 2     -1.31 0.045751634
#> 3     -1.11 0.037037037
#> 4     -0.91 0.038126362
#> 5     -0.70 0.046477850
#> 6     -0.56 0.036310821
#> 7     -0.40 0.032316630
#> 8     -0.28 0.041757444
#> 9     -0.21 0.038489470
#> 10    -0.07 0.053740015
#> 11     0.04 0.210239651
#> 12     0.16 0.039215686
#> 13     0.28 0.038852578
#> 14     0.40 0.037400145
#> 15     0.51 0.055192447
#> 16     0.67 0.042846768
#> 17     0.84 0.053013798
#> 18     1.04 0.057371097
#> 19     1.29 0.049019608
#> 20     1.62 0.002904866
splitpoints(fit, var = "mod_scale", round = 2)
#>    location  proportion
#> 1      0.05 0.067164179
#> 2      0.09 0.045895522
#> 3      0.14 0.042910448
#> 4      0.19 0.043656716
#> 5      0.24 0.106343284
#> 6      0.29 0.041044776
#> 7      0.33 0.041044776
#> 8      0.37 0.053358209
#> 9      0.43 0.045895522
#> 10     0.48 0.042537313
#> 11     0.54 0.069776119
#> 12     0.59 0.101865672
#> 13     0.63 0.044776119
#> 14     0.67 0.042164179
#> 15     0.72 0.038432836
#> 16     0.78 0.039179104
#> 17     0.82 0.036940299
#> 18     0.86 0.038805970
#> 19     0.90 0.052611940
#> 20     0.95 0.005597015
# }
```
