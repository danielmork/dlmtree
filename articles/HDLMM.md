# HDLMM

This vignette demonstrates the implementation of heterogeneous treed
distributed lag mixture model (HDLMM).

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

# Gaussian
hdlmm.fit <- dlmtree(formula = bwgaz ~ ChildSex + MomAge + MomPriorBMI +
                      Race + Hispanic + SmkAny + EstMonthConcept,
                    data = sbd_cov,
                    exposure.data = sbd_exp,
                    family = "gaussian",
                    dlm.type = "linear",
                    mixture = TRUE,
                    het = TRUE,
                    control.het = list(
                      modifiers = c("ChildSex", "MomAge", "MomPriorBMI", "SmkAny"),
                      modifier.splits = 10),
                    control.mcmc = list(n.burn = 2500, n.iter = 10000, n.thin = 5))
```

    #> Preparing data...
    #> 
    #> Running HDLMM:
    #> Burn-in % complete 
    #> [0--------25--------50--------75--------100]
    #>  ''''''''''''''''''''''''''''''''''''''''''
    #> MCMC iterations (est time: 15 minutes)
    #> [0--------25--------50--------75--------100]
    #>  ''''''''''''''''''''''''''''''''''''''''''
    #> Compiling results...

### Model fit summary

``` r

hdlmm.sum <- summary(hdlmm.fit)
print(hdlmm.sum)
```

    #> ---
    #> HDLMM summary
    #> 
    #> Model run info:
    #> - bwgaz ~ ChildSex + MomAge + MomPriorBMI + Race + Hispanic + SmkAny + EstMonthConcept 
    #> - family: gaussian 
    #> - 20 trees
    #> - 2500 burn-in iterations
    #> - 10000 post-burn iterations
    #> - 5 thinning factor
    #> - 5 exposures measured at 37 time points
    #> - 10 two-way interactions (no-self interactions)
    #> - 0.5 modifier sparsity prior
    #> - 1 exposure sparsity prior
    #> - 0.95 confidence level
    #> 
    #> Fixed effects:
    #>                        Mean  Lower  Upper
    #> *(Intercept)          1.515  0.996  1.974
    #>  ChildSexM           -0.454 -1.114  0.274
    #>  MomAge               0.000 -0.003  0.004
    #> *MomPriorBMI         -0.021 -0.025 -0.017
    #>  RaceAsianPI          0.026 -0.094  0.139
    #>  RaceBlack            0.035 -0.092  0.157
    #>  Racewhite            0.016 -0.100  0.129
    #> *HispanicNonHispanic  0.255  0.232  0.277
    #> *SmkAnyY             -0.387 -0.444 -0.251
    #> *EstMonthConcept2     0.118  0.049  0.189
    #> *EstMonthConcept3     0.216  0.113  0.318
    #> *EstMonthConcept4     0.311  0.177  0.443
    #> *EstMonthConcept5     0.421  0.268  0.574
    #> *EstMonthConcept6     0.409  0.246  0.577
    #> *EstMonthConcept7     0.446  0.279  0.608
    #> *EstMonthConcept8     0.427  0.277  0.579
    #> *EstMonthConcept9     0.482  0.355  0.606
    #> *EstMonthConcept10    0.340  0.234  0.447
    #> *EstMonthConcept11    0.223  0.138  0.309
    #>  EstMonthConcept12    0.049 -0.014  0.109
    #> ---
    #> * = CI does not contain zero
    #> 
    #> Modifiers:
    #>                PIP
    #> ChildSex    1.0000
    #> MomAge      0.8550
    #> MomPriorBMI 1.0000
    #> SmkAny      0.0905
    #> ---
    #> PIP = Posterior inclusion probability
    #> 
    #> residual standard errors: 0.02
    #> ---
    #> To obtain exposure effect estimates, use the 'shiny(fit)' function.

### Launching Shiny app

``` r

# shiny(hdlmm.fit)
```
