# HDLM

This vignette demonstrates the implementation of heterogeneous treed
distributed lag model (HDLM). More details can be found in Mork et
al. (2024) \<doi:
[10.1080/01621459.2023.2258595](https://doi.org/10.1080/01621459.2023.2258595)\>.

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

hdlm.fit <- dlmtree(formula = bwgaz ~ ChildSex + MomAge + MomPriorBMI +
                      Race + Hispanic + SmkAny + EstMonthConcept,
                    data = sbd_cov,
                    exposure.data = sbd_exp[["PM25"]],
                    family = "gaussian",
                    dlm.type = "linear",
                    het = TRUE,
                    control.het = list(
                      modifiers = c("ChildSex", "MomAge", "MomPriorBMI", "SmkAny"),
                      modifier.splits = 10),
                    control.mcmc = list(n.burn = 2500, n.iter = 10000, n.thin = 5))
#> Preparing data...
#> 
#> Running shared HDLM:
#> Burn-in % complete 
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> MCMC iterations (est time: 4.1 minutes)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
```

### Model fit summary

``` r

hdlm.sum <- summary(hdlm.fit)
print(hdlm.sum)
#> ---
#> HDLM summary
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
#> - 0.5 modifier sparsity prior
#> - 0.95 confidence level
#> 
#> Fixed effects:
#>                        Mean  Lower  Upper
#> *(Intercept)          1.295  0.949  1.620
#>  ChildSexM            0.115 -0.315  0.551
#>  MomAge               0.000 -0.002  0.002
#> *MomPriorBMI         -0.021 -0.022 -0.019
#>  RaceAsianPI          0.045 -0.079  0.172
#>  RaceBlack            0.055 -0.068  0.182
#>  Racewhite            0.035 -0.083  0.155
#> *HispanicNonHispanic  0.255  0.233  0.277
#>  SmkAnyY             -0.373 -0.558  0.165
#>  EstMonthConcept2    -0.051 -0.106  0.005
#> *EstMonthConcept3    -0.138 -0.201 -0.071
#> *EstMonthConcept4    -0.210 -0.274 -0.144
#> *EstMonthConcept5    -0.198 -0.252 -0.142
#> *EstMonthConcept6    -0.201 -0.254 -0.148
#>  EstMonthConcept7    -0.030 -0.083  0.024
#> *EstMonthConcept8     0.151  0.089  0.212
#> *EstMonthConcept9     0.393  0.330  0.459
#> *EstMonthConcept10    0.379  0.315  0.441
#> *EstMonthConcept11    0.332  0.275  0.387
#> *EstMonthConcept12    0.135  0.085  0.183
#> ---
#> * = CI does not contain zero
#> 
#> Modifiers:
#>                PIP
#> ChildSex    1.0000
#> MomAge      0.0645
#> MomPriorBMI 0.0995
#> SmkAny      0.2110
#> ---
#> PIP = Posterior inclusion probability
#> 
#> residual standard errors: 0.009
#> ---
#> To obtain exposure effect estimates, use the 'shiny(fit)' function.
```

### Launching Shiny app

``` r

# shiny(hdlm.fit)
```
