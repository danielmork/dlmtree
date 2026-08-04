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
#> MCMC iterations (est time: 28 seconds)
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
#> (Intercept)           0.163 -0.863  1.194
#> *ChildSexM           -2.106 -2.127 -2.085
#> MomAge                0.001 -0.001  0.002
#> *MomPriorBMI         -0.021 -0.023 -0.019
#> RaceAsianPI           0.025 -0.095  0.153
#> RaceBlack             0.032 -0.098  0.155
#> Racewhite             0.012 -0.107  0.131
#> *HispanicNonHispanic  0.255  0.232  0.278
#> *SmkAnyY             -0.398 -0.446 -0.350
#> *EstMonthConcept2     0.118  0.034  0.200
#> *EstMonthConcept3     0.234  0.104  0.367
#> *EstMonthConcept4     0.371  0.209  0.535
#> *EstMonthConcept5     0.498  0.319  0.677
#> *EstMonthConcept6     0.450  0.277  0.624
#> *EstMonthConcept7     0.386  0.217  0.549
#> *EstMonthConcept8     0.235  0.071  0.397
#> *EstMonthConcept9     0.259  0.103  0.421
#> *EstMonthConcept10    0.152  0.013  0.291
#> *EstMonthConcept11    0.119  0.014  0.227
#> EstMonthConcept12     0.016 -0.061  0.092
#> ---
#> * = CI does not contain zero
#> 
#> DLNM effect:
#> range = [-0.041, 0.063]
#> signal-to-noise = 0.41
#> critical windows: 4-6,10-34 
#> 
#> residual standard errors: 0.008
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
