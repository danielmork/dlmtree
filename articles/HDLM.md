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
#> MCMC iterations (est time: 3.4 minutes)
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
#> *(Intercept)          1.269  0.922  1.623
#>  ChildSexM            0.124 -0.312  0.552
#>  MomAge               0.001 -0.003  0.005
#> *MomPriorBMI         -0.021 -0.023 -0.018
#>  RaceAsianPI          0.047 -0.076  0.169
#>  RaceBlack            0.055 -0.070  0.180
#>  Racewhite            0.036 -0.090  0.150
#> *HispanicNonHispanic  0.255  0.232  0.276
#> *SmkAnyY             -0.406 -0.458 -0.356
#>  EstMonthConcept2    -0.048 -0.101  0.007
#> *EstMonthConcept3    -0.125 -0.185 -0.064
#> *EstMonthConcept4    -0.194 -0.255 -0.133
#> *EstMonthConcept5    -0.192 -0.245 -0.138
#> *EstMonthConcept6    -0.198 -0.248 -0.148
#>  EstMonthConcept7    -0.035 -0.087  0.017
#> *EstMonthConcept8     0.149  0.092  0.207
#> *EstMonthConcept9     0.399  0.334  0.464
#> *EstMonthConcept10    0.395  0.333  0.454
#> *EstMonthConcept11    0.347  0.292  0.398
#> *EstMonthConcept12    0.143  0.097  0.193
#> ---
#> * = CI does not contain zero
#> 
#> Modifiers:
#>                PIP
#> ChildSex    1.0000
#> MomAge      0.9485
#> MomPriorBMI 0.9225
#> SmkAny      0.2005
#> ---
#> PIP = Posterior inclusion probability
#> 
#> residual standard errors: 0.004
#> ---
#> To obtain exposure effect estimates, use the 'shiny(fit)' function.
```

### Launching Shiny app

``` r

# shiny(hdlm.fit)
```
