# Global variables for plots and 'shiny' app
globalVariables(c("lower", "upper", "est", "Modifier", "PIP", "location", "proportion", 
                  "est.lower", "est.upper", "Var1", "Week", "Obs", "Var2",
                  "%>%", ":=", "Iter", "Rule", "Tree", "bind_rows", "coExp", "exposureCov", 
                  "group_by", "mutate", "pivot_longer", "pm25Exposures", "zinbCo", "pull", "sample_n", 
                  "str_detect", "summarize", "X", "CIMin", "CIMax", "Est", "x", "y", "Effect",
                  "CW", "Tmin", "Tmax", "Xmin", "Xmax", "SD", "vals", "lower", "upper", "PredVal",
                  "fit", "iteration", "effect", "value", "acf", "reshape", "tree", "size", "success",
                  "decision", "Row.num", "Count", "Exposure", "V1", "V2"))

#' print
#' 
#' @description print generic function for S3method
#' 
#' @param x An object of class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone', 
#' representing a fitted model using dlmtree(); or a summary object produced by applying 
#' summary() to one of these model objects.
#' @param digits number of decimal places to round the numeric values to
#' @param ... additional parameters
#' 
#' @return
#' For a fitted model object, prints an assorted model output including model formula call and available methods.
#' For a summary object, prints a summary output of a model fit in the R console.
#' @export print
print <- function(x, ...){
  UseMethod("print")
}


#' summary
#' 
#' @description summary generic function for S3method
#' 
#' @param x an object of class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone'
#' @param conf.level confidence level for computation of credible intervals
#' @param marginalize value(s) for calculating marginal DLMs, defaults to "mean", 
#' can also specify a percentile from 1-99 for all other exposures, or
#' a named vector with specific values for each exposure
#' @param log10BF.crit Bayes Factor criteria for selecting exposures and interactions, 
#' such that log10(BayesFactor) > x. Default = 0.5.
#' @param pred.at numerical vector of exposure values to make predictions for 
#' at each time period
#' @param cenval scalar exposure value that acts as a reference point for 
#' predictions at all other exposure values
#' @param exposure.se scalar smoothing factor, if different from model
#' @param mcmc keep all mcmc iterations (large memory requirement)
#' @param verbose show progress in console
#' @param ... additional parameters
#' 
#' @returns list of summary outputs of the model fit
#' 
#' @export summary
summary <- function(x, conf.level = 0.95, ...){
  UseMethod("summary")
}


#' predict
#' 
#' @description predict generic function for S3method
#' 
#' @param x fitted dlmtree model with class 'hdlm', 'hdlmm'
#' @param new.data new data frame which contains the same covariates and modifiers used to fit the model
#' @param new.exposure.data new data frame/list which contains the same length of exposure lags used to fit the model
#' @param ci.level credible interval level for posterior predictive distribution
#' @param type type of prediction: "response" (default) or "waic". "waic" must be specified with `outcome` parameter
#' @param outcome outcome required for WAIC calculation
#' @param fixed.idx fixed index
#' @param est.dlm flag for estimating dlm effect
#' @param verbose TRUE (default) or FALSE: print output
#' @param ... not used
#'
#' @returns list with the following elements:
#' \describe{
#'  \item{ztg}{posterior predictive mean of fixed effect}
#'  \item{ztg.lims}{lower/upper bound of posterior predictive distribution of fixed effect}
#'  \item{dlmest}{estimated exposure effect}
#'  \item{dlmest.lower}{lower bound of estimated exposure effect}
#'  \item{dlmest.upper}{upper bound of estimated exposure effect}
#'  \item{fhat}{posterior predictive mean of exposure effect}
#'  \item{fhat.lims}{lower/upper bound of posterior predictive distribution of exposure effect}
#'  \item{y}{posterior predictive mean}
#'  \item{y.lims}{lower/upper bound of posterior predictive distribution}
#' }
#' @export predict
predict <- function(x, 
                    new.data, 
                    new.exposure.data,
                    ci.level = 0.95, 
                    type = "response", 
                    outcome = NULL,
                    fixed.idx = list(), 
                    est.dlm = FALSE, 
                    verbose = TRUE,
                    ...){
  UseMethod("predict")
}


#' diagnose
#' 
#' @description diagnose generic function for S3method
#' 
#' @param x a summary object resulting from summary() applied to an object of 
#' class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone'
#' @param ... not used.
#' 
#' @returns shiny interface for assessing model convergence. The interface includes 
#' tabs for MCMC diagnostics such as trace plots, density plots, and convergence measures 
#' for distributed lag effects, DLM tree sizes, and hyperparameters.
#' @export diagnose
diagnose <- function(x, ...){
  UseMethod("diagnose")
}


#' shiny
#'
#' @description shiny generic function for S3method
#'
#' @param fit object of class 'hdlm', 'hdlmm' to which S3method is applied
#'
#' @returns shiny interface for further analysis on heterogeneous analyses. 
#' The interface includes tabs for modifier selection, personalized exposure 
#' effects and subgroup-specific effects.
#' @export shiny
shiny <- function(fit) {
  UseMethod("shiny")
}