# TDLNM

This vignette demonstrates the implementation of treed distributed lag
non-linear model (TDLNM). More details can be found in Mork and Wilson
(2021) \<doi:
[10.1093/biostatistics/kxaa051](https://doi.org/10.1093/biostatistics/kxaa051)\>.

``` r

library(dlmtree)
library(dplyr)
set.seed(1)
```

### Load data

Simulated data is available on
[GitHub](https://github.com/danielmork/dlmtree/tree/master/vignettes/articles).
It can be loaded with the following code.

``` r

sbd_dlmtree <- get_sbd_dlmtree()
```

### Data preparation

``` r

# Response and covariates
sbd_cov <- sbd_dlmtree %>% 
            select(bwgaz, ChildSex, MomAge, GestAge, MomPriorBMI, Race,
                    Hispanic, MomEdu, SmkAny, Marital, Income,
                    EstDateConcept, EstMonthConcept, EstYearConcept)

# Exposure data
sbd_exp <- list(PM25 = sbd_dlmtree %>% select(starts_with("pm25_")),
                TEMP = sbd_dlmtree %>% select(starts_with("temp_")),
                SO2 = sbd_dlmtree %>% select(starts_with("so2_")),
                CO = sbd_dlmtree %>% select(starts_with("co_")),
                NO2 = sbd_dlmtree %>% select(starts_with("no2_")))
sbd_exp <- sbd_exp %>% lapply(as.matrix)
```

### Fitting the model

``` r

tdlnm.fit <- dlmtree(formula = bwgaz ~ ChildSex + MomAge + MomPriorBMI + 
                       Race + Hispanic + SmkAny + EstMonthConcept,
                     data = sbd_cov,
                     exposure.data = sbd_exp[["TEMP"]],
                     dlm.type = "nonlinear",
                     family = "gaussian",
                     control.tdlnm = list(exposure.splits = 20),
                     control.mcmc = list(n.burn = 2500, n.iter = 10000, n.thin = 5))
#> Preparing data...
#> 
#> Running TDLNM:
#> Burn-in % complete 
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> MCMC iterations (est time: 16 seconds)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
```

### Model fit summary

``` r

tdlnm.sum <- summary(tdlnm.fit)
#> Centered DLNM at exposure value 0
print(tdlnm.sum)
#> ---
#> TDLNM summary
#> 
#> Model run info:
#> - bwgaz ~ ChildSex + MomAge + MomPriorBMI + Race + Hispanic + SmkAny + EstMonthConcept 
#> - sample size: 10,000 
#> - family: gaussian 
#> - 20 trees
#> - 2500 burn-in iterations
#> - 10000 post-burn iterations
#> - 5 thinning factor
#> - exposure measured at 37 time points
#> - 0.95 confidence level
#> 
#> Fixed effect coefficients:
#>                        Mean  Lower  Upper
#> (Intercept)           0.186 -2.892  3.193
#> *ChildSexM           -2.106 -2.127 -2.085
#> MomAge                0.001 -0.001  0.002
#> *MomPriorBMI         -0.021 -0.023 -0.019
#> RaceAsianPI           0.026 -0.101  0.156
#> RaceBlack             0.030 -0.102  0.160
#> Racewhite             0.012 -0.113  0.138
#> *HispanicNonHispanic  0.256  0.232  0.278
#> *SmkAnyY             -0.397 -0.442 -0.351
#> *EstMonthConcept2     0.123  0.036  0.203
#> *EstMonthConcept3     0.238  0.105  0.362
#> *EstMonthConcept4     0.371  0.211  0.524
#> *EstMonthConcept5     0.500  0.329  0.676
#> *EstMonthConcept6     0.453  0.276  0.630
#> *EstMonthConcept7     0.391  0.220  0.567
#> *EstMonthConcept8     0.241  0.074  0.413
#> *EstMonthConcept9     0.268  0.106  0.436
#> *EstMonthConcept10    0.160  0.017  0.297
#> *EstMonthConcept11    0.123  0.017  0.230
#> EstMonthConcept12     0.014 -0.062  0.086
#> ---
#> * = CI does not contain zero
#> 
#> DLNM effect:
#> range = [-0.047, 0.067]
#> signal-to-noise = 0.411
#> critical windows: 1-7,11-34 
#> 
#> residual standard errors: 0.011
```

### Exposure-time surface

``` r

plot(tdlnm.sum, 
     main = "Plot title", 
     xlab = "Time axis label", 
     ylab = "Exposure-concentration axis label", 
     flab = "Effect color label")
```

![](TDLNM_files/figure-html/tdlnm.plot-1.png)

### Slicing on exposure-concentration

``` r

# slicing on exposure-concentration
plot(tdlnm.sum, plot.type = "slice", val = 1, main = "Slice at concentration 1") 
```

![](TDLNM_files/figure-html/unnamed-chunk-3-1.png)

``` r

plot(tdlnm.sum, plot.type = "slice", val = 2, main = "Slice at concentration 2")
```

![](TDLNM_files/figure-html/unnamed-chunk-3-2.png)

### Slicing on time lag

``` r

# slicing on exposure-concentration
plot(tdlnm.sum, plot.type = "slice", time = 7, main = "Slice at time 7")
```

![](TDLNM_files/figure-html/unnamed-chunk-4-1.png)

``` r

plot(tdlnm.sum, plot.type = "slice", time = 15, main = "Slice at time 15")
```

![](TDLNM_files/figure-html/unnamed-chunk-4-2.png)

``` r

plot(tdlnm.sum, plot.type = "slice", time = 33, main = "Slice at time 33")
```

![](TDLNM_files/figure-html/unnamed-chunk-4-3.png)

### different plot.type options

``` r

# Standard error, credible intervals
plot(tdlnm.sum, plot.type = "se", main = "Standard error")  
```

![](TDLNM_files/figure-html/unnamed-chunk-5-1.png)

``` r

plot(tdlnm.sum, plot.type = "ci-min", main = "Credible interval lower bound")
```

![](TDLNM_files/figure-html/unnamed-chunk-5-2.png)

``` r

plot(tdlnm.sum, plot.type = "ci-max", main = "Credible interval upper bound")
```

![](TDLNM_files/figure-html/unnamed-chunk-5-3.png)

``` r

# Cumulative effect and significance
plot(tdlnm.sum, plot.type = "cumulative", main = "Cumulative effect per exposure-concentration")
```

![](TDLNM_files/figure-html/unnamed-chunk-6-1.png)

``` r

plot(tdlnm.sum, plot.type = "effect", main = "Significant effects with directions")
```

![](TDLNM_files/figure-html/unnamed-chunk-6-2.png)
