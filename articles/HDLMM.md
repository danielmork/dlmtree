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
    #> MCMC iterations (est time: 17 minutes)
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
    #> *(Intercept)          1.537  1.063  1.961
    #>  ChildSexM           -0.458 -0.988  0.042
    #>  MomAge               0.000 -0.002  0.004
    #> *MomPriorBMI         -0.021 -0.025 -0.017
    #>  RaceAsianPI          0.025 -0.104  0.155
    #>  RaceBlack            0.034 -0.098  0.160
    #>  Racewhite            0.015 -0.108  0.138
    #> *HispanicNonHispanic  0.255  0.233  0.278
    #> *SmkAnyY             -0.382 -0.443 -0.154
    #> *EstMonthConcept2     0.120  0.052  0.193
    #> *EstMonthConcept3     0.219  0.118  0.316
    #> *EstMonthConcept4     0.314  0.181  0.445
    #> *EstMonthConcept5     0.425  0.276  0.575
    #> *EstMonthConcept6     0.413  0.256  0.574
    #> *EstMonthConcept7     0.450  0.298  0.608
    #> *EstMonthConcept8     0.427  0.288  0.572
    #> *EstMonthConcept9     0.478  0.351  0.608
    #> *EstMonthConcept10    0.338  0.225  0.452
    #> *EstMonthConcept11    0.224  0.136  0.313
    #>  EstMonthConcept12    0.050 -0.008  0.111
    #> ---
    #> * = CI does not contain zero
    #> 
    #> Modifiers:
    #>                PIP
    #> ChildSex    1.0000
    #> MomAge      0.8865
    #> MomPriorBMI 1.0000
    #> SmkAny      0.1630
    #> ---
    #> PIP = Posterior inclusion probability
    #> 
    #> residual standard errors: 0.008
    #> ---
    #> To obtain exposure effect estimates, use the 'shiny(fit)' function.

### Launching Shiny app

``` r

# shiny(hdlmm.fit)
```
