# TDLM

This vignette demonstrates the implementation of treed distributed lag
model (TDLM). More details can be found in Mork and Wilson (2023) \<doi:
[10.1111/biom.13568](https://doi.org/10.1111/biom.13568)\>.

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

tdlm.fit <- dlmtree(formula = bwgaz ~ ChildSex + MomAge + MomPriorBMI +
                      Race + Hispanic + SmkAny + EstMonthConcept,
                    data = sbd_cov,
                    exposure.data = sbd_exp[["PM25"]], # A single numeric matrix
                    family = "gaussian",
                    dlm.type = "linear",
                    control.mcmc = list(n.burn = 2500, n.iter = 10000, n.thin = 5))
```

    #> Preparing data...
    #> 
    #> Running TDLM:
    #> Burn-in % complete 
    #> [0--------25--------50--------75--------100]
    #>  ''''''''''''''''''''''''''''''''''''''''''
    #> MCMC iterations (est time: 16 seconds)
    #> [0--------25--------50--------75--------100]
    #>  ''''''''''''''''''''''''''''''''''''''''''
    #> Compiling results...

### Model fit summary

``` r

tdlm.sum <- summary(tdlm.fit)
print(tdlm.sum)
```

    #> ---
    #> TDLM summary
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
    #> *(Intercept)          2.409  2.135  2.679
    #> *ChildSexM           -2.106 -2.126 -2.086
    #> MomAge                0.000 -0.001  0.002
    #> *MomPriorBMI         -0.021 -0.022 -0.019
    #> RaceAsianPI           0.044 -0.082  0.173
    #> RaceBlack             0.055 -0.077  0.186
    #> Racewhite             0.035 -0.084  0.160
    #> *HispanicNonHispanic  0.254  0.231  0.276
    #> *SmkAnyY             -0.404 -0.451 -0.358
    #> EstMonthConcept2     -0.053 -0.113  0.011
    #> *EstMonthConcept3    -0.149 -0.215 -0.084
    #> *EstMonthConcept4    -0.232 -0.299 -0.162
    #> *EstMonthConcept5    -0.206 -0.265 -0.150
    #> *EstMonthConcept6    -0.205 -0.259 -0.149
    #> EstMonthConcept7     -0.034 -0.091  0.023
    #> *EstMonthConcept8     0.140  0.077  0.205
    #> *EstMonthConcept9     0.387  0.317  0.454
    #> *EstMonthConcept10    0.367  0.300  0.431
    #> *EstMonthConcept11    0.329  0.271  0.387
    #> *EstMonthConcept12    0.127  0.075  0.181
    #> ---
    #> * = CI does not contain zero
    #> 
    #> DLM effect:
    #> range = [-0.019, 0.008]
    #> signal-to-noise = 0.022
    #> critical windows: 11-20 
    #>              Mean  Lower  Upper
    #> Period 1    0.002 -0.005  0.015
    #> Period 2    0.000 -0.007  0.009
    #> Period 3   -0.002 -0.010  0.004
    #> Period 4   -0.003 -0.011  0.002
    #> Period 5   -0.002 -0.008  0.004
    #> Period 6   -0.001 -0.006  0.005
    #> Period 7   -0.001 -0.006  0.005
    #> Period 8   -0.002 -0.008  0.004
    #> Period 9   -0.003 -0.013  0.004
    #> Period 10  -0.003 -0.014  0.005
    #> *Period 11 -0.016 -0.024 -0.007
    #> *Period 12 -0.017 -0.024 -0.010
    #> *Period 13 -0.018 -0.024 -0.012
    #> *Period 14 -0.017 -0.022 -0.011
    #> *Period 15 -0.017 -0.022 -0.011
    #> *Period 16 -0.017 -0.022 -0.011
    #> *Period 17 -0.018 -0.023 -0.012
    #> *Period 18 -0.019 -0.030 -0.013
    #> *Period 19 -0.019 -0.028 -0.011
    #> *Period 20 -0.015 -0.024 -0.002
    #> Period 21  -0.007 -0.019  0.002
    #> Period 22  -0.003 -0.011  0.005
    #> Period 23  -0.003 -0.012  0.002
    #> Period 24  -0.002 -0.009  0.004
    #> Period 25  -0.001 -0.006  0.005
    #> Period 26  -0.001 -0.005  0.005
    #> Period 27  -0.001 -0.006  0.004
    #> Period 28  -0.001 -0.007  0.003
    #> Period 29  -0.002 -0.008  0.003
    #> Period 30  -0.002 -0.009  0.003
    #> Period 31  -0.002 -0.009  0.004
    #> Period 32  -0.001 -0.008  0.006
    #> Period 33   0.002 -0.005  0.009
    #> Period 34   0.003 -0.004  0.011
    #> Period 35   0.005 -0.002  0.013
    #> Period 36   0.007  0.000  0.016
    #> Period 37   0.008  0.000  0.017
    #> ---
    #> * = CI does not contain zero
    #> 
    #> residual standard errors: 0.01
    #> ---

### Exposure effect

``` r

plot(tdlm.sum, 
     main = "Plot title", 
     xlab = "Time axis label", 
     ylab = "Exposure effect axis label")
```

![](TDLM_files/figure-html/tdlm.plot-1.png)
