#' splitpoints
#'
#' @title Determines split points for continuous modifiers
#' @description Method for determining split points for continuous modifiers
#'
#' @param object An object of class dlmtree with DLM type hdlm & hdlmm
#' @param var The name of a continuous variable for which the split points will be reported
#' @param round The number of decimal places to round the variable (var) to. No rounding occurs if round=NULL (default)
#' For positive integer values of round, the variable will be rounded and split points will be reported at the resulting level
#'
#' @examples
#' # Split points with HDLM 
#' D <- sim.hdlmm(sim = "B", n = 1000)
#' fit <- dlmtree(y ~ ., 
#'                data = D$dat,
#'                exposure.data = D$exposures,
#'                dlm.type = "linear",
#'                family = "gaussian",
#'                het = TRUE)
#' splitpoints(fit, var = "mod_num", round = 2)
#' splitpoints(fit, var = "mod_scale", round = 2)
#'
#' @returns A data frame with split points and the probability that a split point was >= that split point value
#' @export
splitpoints <- function(object, var, round = NULL) 
{
  sp2         <- object$TreeStructs$Rule[!duplicated(object$TreeStructs[,2:4])]
  treeRules   <- object$TreeStructs %>% group_by(Iter, Tree) %>% 
    summarize(Rules = paste0(Rule, collapse = " & "))
  splitRules2 <- table(do.call(c, lapply(strsplit(treeRules$Rules, " & ", T), unique)))
  
  # check if continuous
  categorical <- length(splitRules2[grepl(var, names(splitRules2)) & 
                        grepl("%in%", names(splitRules2))])>0
                        
  if (categorical) {
    stop("var is categorical. Split points only works with continuous variables")
  } else {
    varsp <- sort(splitRules2[grepl(var, names(splitRules2)) & 
                                grepl(">=", names(splitRules2))]) /
              sum(sort(splitRules2[grepl(var, names(splitRules2)) & 
                                    grepl(">=", names(splitRules2))]))
    
    splits <- data.frame(loc = as.numeric(sapply(strsplit(names(varsp), " >= ", T), function(i) i[2])),
                         val = as.numeric(varsp))
    
    colnames(splits) <- c("location","proportion")
    
    if (nrow(splits) == 0) {
      cat("There are either no splits")
    }
    
    if (!is.null(round)) {
      if (is.numeric(splits$location)) {
        splits$location <- round(splits$location,round)
      } else {
        cat("rounding is not permitted when var is not numeric.")
      }
    }
    
    if (is.numeric(splits$location)) {
      splits <- aggregate(proportion ~ location, data=splits, sum)
    }
    
    return(splits)
  }
  
}

