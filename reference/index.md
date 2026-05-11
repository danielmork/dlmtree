# Package index

## Model fitting

- [`dlmtree()`](https://danielmork.github.io/dlmtree/reference/dlmtree.md)
  : Fit tree structured distributed lag models
- [`dlmtree.control.diagnose()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.diagnose.md)
  : Diagnostic control settings for dlmtree model fitting
- [`dlmtree.control.family()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.family.md)
  : Family control settings for dlmtree model fitting
- [`dlmtree.control.het()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.het.md)
  : Control settings for dlmtree model fitting, when used for
  heterogeneous models
- [`dlmtree.control.hyper()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.hyper.md)
  : Hyperparameter control settings for dlmtree model fitting
- [`dlmtree.control.mcmc()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.mcmc.md)
  : MCMC control settings for dlmtree model fitting
- [`dlmtree.control.mix()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.mix.md)
  : Control settings for dlmtree model fitting, when used for mixture
  models
- [`dlmtree.control.monotone()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.monotone.md)
  : Control settings for dlmtree model fitting, when used for monotone
  model
- [`dlmtree.control.tdlnm()`](https://danielmork.github.io/dlmtree/reference/dlmtree.control.tdlnm.md)
  : Control settings for dlmtree model fitting, when used for TDLNM
- [`diagnose()`](https://danielmork.github.io/dlmtree/reference/diagnose.md)
  : diagnose

## Model summary

Obtain results from the model fit

- [`combine.models()`](https://danielmork.github.io/dlmtree/reference/combine.models.md)
  : Combines information from DLMs of single exposure
- [`combine.models.tdlmm()`](https://danielmork.github.io/dlmtree/reference/combine.models.tdlmm.md)
  : Combines information from DLMs of mixture exposures.
- [`summary()`](https://danielmork.github.io/dlmtree/reference/summary.md)
  : summary
- [`estDLM()`](https://danielmork.github.io/dlmtree/reference/estDLM.md)
  : Calculates subgroup-specific lag effects for heterogeneous models
- [`ppRange()`](https://danielmork.github.io/dlmtree/reference/ppRange.md)
  : Makes a 'pretty' output of a group of numbers
- [`predict()`](https://danielmork.github.io/dlmtree/reference/predict.md)
  : predict
- [`adj_coexposure()`](https://danielmork.github.io/dlmtree/reference/adj_coexposure.md)
  : Adjusting for expected changes in co-exposure (TDLMM)
- [`print()`](https://danielmork.github.io/dlmtree/reference/print.md) :
  print

## Visualization

Visualize distributed lag functions

- [`plot(`*`<summary.monotone>`*`)`](https://danielmork.github.io/dlmtree/reference/plot.summary.monotone.md)
  : Returns variety of plots for model summary of class 'monotone'
- [`plot(`*`<summary.tdlm>`*`)`](https://danielmork.github.io/dlmtree/reference/plot.summary.tdlm.md)
  : Plots a distributed lag function for model summary of 'tdlm'
- [`plot(`*`<summary.tdlmm>`*`)`](https://danielmork.github.io/dlmtree/reference/plot.summary.tdlmm.md)
  : Plots DLMMs for model summary of class 'tdlmm'
- [`plot(`*`<summary.tdlnm>`*`)`](https://danielmork.github.io/dlmtree/reference/plot.summary.tdlnm.md)
  : Returns variety of plots for model summary of class 'tdlnm'

## HDLM & HDLMM

Useful functions for HDLM & HDLMM

- [`pip()`](https://danielmork.github.io/dlmtree/reference/pip.md) :
  Calculates posterior inclusion probabilities (PIPs) for modifiers in
  HDLM & HDLMM
- [`shiny()`](https://danielmork.github.io/dlmtree/reference/shiny.md) :
  shiny
- [`splitpoints()`](https://danielmork.github.io/dlmtree/reference/splitpoints.md)
  : Determines split points for continuous modifiers

## Simulation

Simulate datasets

- [`sim.hdlmm()`](https://danielmork.github.io/dlmtree/reference/sim.hdlmm.md)
  : Creates simulated data for HDLM & HDLMM
- [`sim.tdlmm()`](https://danielmork.github.io/dlmtree/reference/sim.tdlmm.md)
  : Creates simulated data for TDLM & TDLMM
- [`sim.tdlnm()`](https://danielmork.github.io/dlmtree/reference/sim.tdlnm.md)
  : Creates simulated data for TDLNM

## Data

Built-in datasets

- [`coExp`](https://danielmork.github.io/dlmtree/reference/coExp.md) :
  Randomly sampled exposure from Colorado counties
- [`exposureCov`](https://danielmork.github.io/dlmtree/reference/exposureCov.md)
  : Exposure covariance structure
- [`pm25Exposures`](https://danielmork.github.io/dlmtree/reference/pm25Exposures.md)
  : PM2.5 Exposure data
- [`zinbCo`](https://danielmork.github.io/dlmtree/reference/zinbCo.md) :
  Time-series exposure data for ZINB simulated data
- [`get_sbd_dlmtree()`](https://danielmork.github.io/dlmtree/reference/get_sbd_dlmtree.md)
  : Download simulated data for dlmtree articles

## Cpp source code & Misc.

- [`cppIntersection()`](https://danielmork.github.io/dlmtree/reference/cppIntersection.md)
  : fast set intersection tool assumes sorted vectors A and B
- [`dlmEst()`](https://danielmork.github.io/dlmtree/reference/dlmEst.md)
  : Calculates the distributed lag effect with DLM matrix for linear
  models.
- [`dlmtreeGPFixedGaussian()`](https://danielmork.github.io/dlmtree/reference/dlmtreeGPFixedGaussian.md)
  : dlmtree model with fixed Gaussian process approach
- [`dlmtreeGPGaussian()`](https://danielmork.github.io/dlmtree/reference/dlmtreeGPGaussian.md)
  : dlmtree model with Gaussian process approach
- [`dlmtreeHDLMGaussian()`](https://danielmork.github.io/dlmtree/reference/dlmtreeHDLMGaussian.md)
  : dlmtree model with shared HDLM approach
- [`dlmtreeHDLMMGaussian()`](https://danielmork.github.io/dlmtree/reference/dlmtreeHDLMMGaussian.md)
  : dlmtree model with HDLMM approach
- [`dlmtreeTDLMFixedGaussian()`](https://danielmork.github.io/dlmtree/reference/dlmtreeTDLMFixedGaussian.md)
  : dlmtree model with fixed Gaussian approach
- [`dlmtreeTDLMNestedGaussian()`](https://danielmork.github.io/dlmtree/reference/dlmtreeTDLMNestedGaussian.md)
  : dlmtree model with nested Gaussian approach
- [`dlmtreeTDLM_cpp()`](https://danielmork.github.io/dlmtree/reference/dlmtreeTDLM_cpp.md)
  : dlmtree model with nested HDLM approach
- [`dlnmEst()`](https://danielmork.github.io/dlmtree/reference/dlnmEst.md)
  : Calculates the distributed lag effect with DLM matrix for non-linear
  models.
- [`dlnmPLEst()`](https://danielmork.github.io/dlmtree/reference/dlnmPLEst.md)
  : Calculates the distributed lag effect with DLM matrix for non-linear
  models.
- [`drawTree()`](https://danielmork.github.io/dlmtree/reference/drawTree.md)
  : Draws a new tree structure
- [`mixEst()`](https://danielmork.github.io/dlmtree/reference/mixEst.md)
  : Calculates the lagged interaction effects with MIX matrix for linear
  models.
- [`monotdlnm_Cpp()`](https://danielmork.github.io/dlmtree/reference/monotdlnm_Cpp.md)
  : dlmtree model with monotone tdlnm approach
- [`ppRange()`](https://danielmork.github.io/dlmtree/reference/ppRange.md)
  : Makes a 'pretty' output of a group of numbers
- [`rcpp_pgdraw()`](https://danielmork.github.io/dlmtree/reference/rcpp_pgdraw.md)
  : Multiple draw polya gamma latent variable for var c\[i\] with size
  b\[i\]
- [`rtmvnorm()`](https://danielmork.github.io/dlmtree/reference/rtmvnorm.md)
  : Truncated multivariate normal sampler, mean mu, cov sigma, truncated
  (0, Inf)
- [`ruleIdx()`](https://danielmork.github.io/dlmtree/reference/ruleIdx.md)
  : Calculates the weights for each modifier rule
- [`scaleModelMatrix()`](https://danielmork.github.io/dlmtree/reference/scaleModelMatrix.md)
  : Centers and scales a matrix
- [`splitPIP()`](https://danielmork.github.io/dlmtree/reference/splitPIP.md)
  : Calculates the posterior inclusion probability (PIP).
- [`tdlmm_Cpp()`](https://danielmork.github.io/dlmtree/reference/tdlmm_Cpp.md)
  : dlmtree model with tdlmm approach
- [`tdlnm_Cpp()`](https://danielmork.github.io/dlmtree/reference/tdlnm_Cpp.md)
  : dlmtree model with tdlnm approach
- [`zeroToInfNormCDF()`](https://danielmork.github.io/dlmtree/reference/zeroToInfNormCDF.md)
  : Integrates (0,inf) over multivariate normal
