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
#> MCMC iterations (est time: 2 seconds)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
splitpoints(fit, var = "mod_num", round = 2)
#>    location   proportion
#> 1     -1.70 0.0432480141
#> 2     -1.32 0.0414827891
#> 3     -1.09 0.0547219771
#> 4     -0.90 0.0247131509
#> 5     -0.73 0.0370697264
#> 6     -0.61 0.0308914387
#> 7     -0.46 0.0273609885
#> 8     -0.34 0.0423654016
#> 9     -0.23 0.0538393645
#> 10    -0.11 0.0661959400
#> 11    -0.02 0.2241835834
#> 12     0.09 0.0723742277
#> 13     0.23 0.0379523389
#> 14     0.35 0.0370697264
#> 15     0.46 0.0423654016
#> 16     0.59 0.0423654016
#> 17     0.79 0.0335392763
#> 18     0.98 0.0494263019
#> 19     1.20 0.0379523389
#> 20     1.63 0.0008826125
splitpoints(fit, var = "mod_scale", round = 2)
#>    location proportion
#> 1      0.05 0.05546751
#> 2      0.10 0.06275753
#> 3      0.15 0.04405705
#> 4      0.19 0.05071315
#> 5      0.23 0.10966719
#> 6      0.28 0.04754358
#> 7      0.32 0.03391442
#> 8      0.37 0.03106181
#> 9      0.42 0.03835182
#> 10     0.47 0.04278922
#> 11     0.51 0.04659271
#> 12     0.55 0.03676704
#> 13     0.61 0.04690967
#> 14     0.66 0.04437401
#> 15     0.71 0.03993661
#> 16     0.76 0.03866878
#> 17     0.81 0.04817750
#> 18     0.85 0.05261490
#> 19     0.91 0.11885895
#> 20     0.95 0.01077655
# }
```
