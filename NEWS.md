# dlmtree 1.1.0
* new feature: `diagnose()` function - launches a Shiny panel for assessing MCMC convergence, which includes:
  * trace plots and density plots for the estimated distributed lag function
  * lag-wise convergence analysis
  * tree update convergence
  * trace plots and density plots for hyperparameters
  * exposure selection (for mixture models)
  * modifier selection (for heterogeneous models)

* update `summary.X()` functions
  * more consistent output formats
  * add options to extract MCMC posterior samples

* documentation improvement
  * add more detailed descriptions in `Value` sections
  * reduced redundancy across help files of S3method

# dlmtree 1.0.1
* change the documentation organization format of the main function using `list()` for input argument
  - lists of arguments are linked to `control.R`
  - remove the prefixes from control arguments
* remove following scripts from the package: `tdlmm.R`, `tdlnm.R`
* update credible interval bands to be transparent

# dlmtree 1.0.0.2
* fix to `shiny.hdlm` and `shiny.hdlmm` to remove bug where no data is selected
* fix to `combine.models.R` to correctly combine tree posterior samples
* fix to `summary.tdlmm` in rowMeans to stop autodrop of dimension

# dlmtree 1.0.0.1
* a split point for Subgroup tab in shiny interface set to the mode of proportion

# dlmtree 1.0.0
* initial CRAN release
* spelling checks using `spell_check_package()` function from `spelling` package
* S3 method `predict.dlmtree` is separately defined as `predict.hdlm` and `predict.hdlmm`
* combined model fitting functions `tdlnm` and `tdlmm` into `dlmtree`

# dlmtree 0.9.0
* Beta prior for monotone CW selection
* Bayes factor selection method
* zero inflated tree method
* monotone DLNM
* precalculation in tdlnm for large speed improvement

# dlmtree 0.8.1.1
* fix to single covariate model
* add functions for posterior analysis with dlmtree class (`pip` and `splitpoints`)

# dlmtree 0.8.1.0
* fix to allow intercept only models
* updated tdlnm algorithm for increased performance (memory restriction for large datasets)

# dlmtree 0.8.0.0
* add HDLM functionality (`dlmtree`)

# dlmtree 0.7.1.1
* fix `summary.tdlmm` if nMix=1 or nExp=1

# dlmtree 0.7.1.0
* combined identity and logistic link into single c++ file
* several code efficiencies
* updated print output of `tdlmm.summary`
* removed prior on dirichlet hyperparameter
* exposure selection based on Beta-Binomial distribution
* `keep.mcmc` in summary

# dlmtree 0.6.2.0
* separated `tdlmm` from `tdlnm` function
* included exposure selection into summary
* prior on dirichlet hyperparameter ~Gamma(1, nTrees/4)
