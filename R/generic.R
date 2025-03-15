# Global variables for plots and 'shiny' app
globalVariables(c("lower", "upper", "est", "Modifier", "PIP", "location", "proportion", 
                  "est.lower", "est.upper", "Var1", "Week", "Obs", "Var2",
                  "%>%", ":=", "Iter", "Rule", "Tree", "bind_rows", "coExp", "exposureCov", 
                  "group_by", "mutate", "pivot_longer", "pm25Exposures", "zinbCo", "pull", "sample_n", 
                  "str_detect", "summarize", "X", "CIMin", "CIMax", "Est", "x", "y", "Effect",
                  "CW", "Tmin", "Tmax", "Xmin", "Xmax", "SD", "vals", "lower", "upper", "PredVal"))

#' print
#' 
#' @description print generic function for S3method
#' 
#' @param x an object of class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone'
#' @param ... not used.
#' 
#' @return assorted model output including model formula call and available methods
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
#' @export summary
summary <- function(x, conf.level = 0.95, ...){
  UseMethod("summary")
}


#' print.summary
#' 
#' @description print.summary generic function for S3method
#' 
#' @param x a summary object resulting from summary() applied to an object of 
#' class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone'
#' @param digits integer number of digits to round
#' @param cw.only print only results for exposures with critical windows
#' @param ... additional parameters
#' 
#' @returns summary output of a model fit printed in the R console
#' @export print.summary
print.summary <- function(x, ...){
  UseMethod("print.summary")
}


#' shiny
#'
#' @description shiny generic function for S3method
#'
#' @param fit an object of class 'hdlm', 'hdlmm' to which S3method is applied
#'
#' @returns 'shiny' interface for further analysis on heterogeneous analyses. 
#' The interface includes tabs for modifier selection, personalized exposure 
#' effects and subgroup-specific effects.
#' @export shiny
shiny <- function(fit) {
  UseMethod("shiny")
}