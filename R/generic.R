# Global variables for plots and shiny
globalVariables(c("lower", "upper", "est", "Modifier", "PIP", "location", "proportion", 
                  "est.lower", "est.upper", "Var1", "Week", "Obs", "Var2",
                  "%>%", ":=", "Iter", "Rule", "Tree", "bind_rows", "coExp", "exposureCov", 
                  "group_by", "mutate", "pivot_longer", "pm25Exposures", "pull", "sample_n", 
                  "str_detect", "summarize", "X", "CIMin", "CIMax", "Est", "x", "y", "Effect", 
                  "CW", "Tmin", "Tmax", "Xmin", "Xmax", "SD", "vals", "lower", "upper", "PredVal"))

#' shiny
#'
#' @description shiny generic function for S3method
#'
#' @param fit an object to which S3method is applied
#'
#' @export shiny
shiny <- function(fit) {
  UseMethod("shiny")
}
