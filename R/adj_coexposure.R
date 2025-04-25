#' adj_coexposure
#'
#' @title Adjusting for expected changes in co-exposure (TDLMM)
#' @description Estimates the marginal effects of an exposure while accounting
#' for expected changes in co-occurring exposures at the same time point.
#' Values of co-occurring exposures are modeled nonlinearly using a spline
#' model with predictions made at the lower an upper values for the exposure
#' of interest.
#' 
#' @param exposure.data Named list of exposure matrices used as input to TDLMM.
#' @param object Model output for TDLMM from dlmtree() function.
#' @param contrast_perc 2-length vector of percentiles or named list 
#' corresponding to lower and upper exposure percentiles of interest.
#' Names must equal list names in 'exposure.data'.
#' @param contrast_exp Named list consisting lower and upper exposure values.
#' This takes precedence over contrast_perc if both inputs are used.
#' @param conf.level Confidence level used for estimating credible intervals.
#' Default is 0.95.
#' @param keep.mcmc If TRUE, return posterior samples.
#' @param verbose TRUE (default) or FALSE: print output
#'
#' @return A list with the following components (or posterior samples if keep.mcmc = TRUE):
#' \item{Name}{vector of exposure names}
#' \item{Time}{integer vector of lags}
#' \item{Effect}{posterior mean of marginal effects}
#' \item{SE}{standard error of the estimate}
#' \item{Lower}{lower bound of credible interval of the marginal effect estimate}
#' \item{Upper}{upper bound of credible interval of the marginal effect estimate}
#' \item{cEffect}{cumulative marginal effects}
#' \item{cLower}{lower bound of credible interval of the cumulative marginal effect}
#' \item{cUpper}{upper bound of credible interval of the cumulative marginal effect}
#' \item{CW}{boolean vector indicating critical window}
#' @importFrom mgcv bam
#' @export

adj_coexposure <- function(exposure.data,
                           object,
                           contrast_perc = c(0.25, 0.75),
                           contrast_exp = list(),
                           conf.level = 0.95,
                           keep.mcmc = FALSE,
                           verbose = TRUE) 
{
  if (!inherits(object, "tdlmm"))
    stop("adj_coexposure is intended to be used with TDLMM")

  ##### Adjusting for changes in co-exposures ######
  exposureDat <- do.call(cbind.data.frame, 
                         lapply(exposure.data, function(e) c(e)))
  expMean <- colMeans(exposureDat) # Empirical means
  if (length(contrast_exp) == 0) {
    if (is.list(contrast_perc)) {
      if (!all.equal(names(contrast_perc), names(exposure.data)))
        stop("named list contrast_perc must have same names as exposure.data")
      contrast_exp <- lapply(names(exposure.data), function(i) quantile(exposure.data[[i]], probs = contrast_perc[[i]]))
      names(contrast_exp) <- names(exposure.data)
    } else {
      contrast_exp <- lapply(exposureDat, function(i) quantile(i, probs = contrast_perc))
    }
  }
    

  if (!all.equal(names(contrast_exp), names(exposure.data)))
    stop("named list contrast_exp must have same names as exposure.data")

  if(verbose){
    print(do.call(cbind.data.frame, contrast_exp))
  }
  
  # predLevels[[predictor]][[response]] =
  #  use level of 'predictor' to estimate 'response'.
  predLevels <- lapply(exposureDat, function(i) list())
  # Loop over pairs of exposures
  for (predExp in names(predLevels)) {
    for (respExp in names(predLevels)) {
      if(verbose){
        cat("\nPredicting", respExp, "at lower/upper values of", predExp)
      }
      if (predExp == respExp) { 
        # If predictor/response exposure same, set to 25/75 percentiles.
        predLevels[[predExp]][[respExp]] <- contrast_exp[[predExp]]
        
      } else { 
        # Predict `response` (respExp) at 25/75 percentiles of `predictor` (predExp).
        # Other methods could be substituted here to achieve prediction
        formula <- as.formula(paste(respExp, "~s(", predExp, ", k = 5, bs = 'ps')"))
        model <- bam(formula, data = exposureDat)
        newdat <- data.frame(contrast_exp[[predExp]])
        names(newdat) <- predExp
        predLevels[[predExp]][[respExp]] <- predict(model, newdata = newdat)
      }
      if(verbose){
        cat(":", predLevels[[predExp]][[respExp]])
      }
    }
  }
  
  
  modelSummary <- summary(object, mcmc = TRUE, verbose = FALSE)
  expDLM <- list() # For output
  
  for (exposure in names(modelSummary$DLM)) {
    # MCMC results for main effect of interest
    mcmc <- modelSummary$DLM[[exposure]]$mcmc * diff(contrast_exp[[exposure]])
    # Add MCMC results for expected changes in co-exposures
    for (coexposure in names(modelSummary$DLM)) {
      if (exposure != coexposure)
        mcmc <- mcmc + 
          modelSummary$DLM[[coexposure]]$mcmc * diff(predLevels[[exposure]][[coexposure]])
    }
    # Add MCMC results due to interactions: 
    # 1) same time interactions, use expected changes in main and co-exposure
    # 2) diff time interactions, use expected change in main, mean in co-exposure
    # 3) same time interactions (b/t co-exp), use expected change in co-exposures
    # 4) diff time interactions (b/t co-exp), cancels (set at mean)
    if (length(modelSummary$MIX) > 0) {
      for (mix in 1:length(modelSummary$MIX)) {
        # Interaction between main and co-exposure
        if (modelSummary$MIX[[mix]]$rows == exposure) {
          mcmc <- mcmc +
            # add same time interactions using expected changes
            t(sapply(1:modelSummary$n.lag, function(k) {
              modelSummary$MIX[[mix]]$mcmc[k,k,]
            })) * 
            (predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][2] * 
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][2] -
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][1] * 
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][1]) +
            # add different time interactions using IQR change in primary exposure
            # (other exposures set to mean)
            t(sapply(1:modelSummary$n.lag, function(k) {
              colSums(modelSummary$MIX[[mix]]$mcmc[k,-k,]) # interaction sum removing time k
            })) * 
            diff(predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]]) *
            expMean[modelSummary$MIX[[mix]]$cols]
          # Interaction between main and co-exposure
        } else if (modelSummary$MIX[[mix]]$cols == exposure) {
          mcmc <- mcmc +
            # add same time interactions using expected changes
            t(sapply(1:modelSummary$n.lag, function(k) {
              modelSummary$MIX[[mix]]$mcmc[k,k,]
            })) * 
            (predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][2] * 
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][2] -
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][1] * 
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][1]) +
            # add different time interactions using IQR change in primary exposure
            # (other exposure set to mean)
            t(sapply(1:modelSummary$n.lag, function(k) {
              colSums(modelSummary$MIX[[mix]]$mcmc[-k,k,]) # interaction sum removing time k
            })) * 
            diff(predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]]) *
            expMean[modelSummary$MIX[[mix]]$rows]
          # Interaction between two co-exposures
        } else {
          mcmc <- mcmc +
            # add same time interactions using expected changes
            t(sapply(1:modelSummary$n.lag, function(k) {
              modelSummary$MIX[[mix]]$mcmc[k,k,]
            })) * 
            (predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][2] * 
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][2] -
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][1] * 
               predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][1])
          # no different time interactions (they are set at mean, so equals zero)
        }
      }
    }
    # Summarize main effects, adjusted for changes in co-exposures
    if (keep.mcmc) {
      expDLM[[exposure]] <- mcmc
    } else {
      expDLM[[exposure]] <- 
        suppressWarnings(
        data.frame("Name" = exposure,
                   "Time" = 1:modelSummary$n.lag,
                   "Effect" = apply(mcmc, 1, mean),
                   "SE" = sqrt(diag(cov(t(mcmc)))),
                   # Credible interval using conf.level
                   "Lower" =  apply(mcmc, 1, quantile, probs = (1 - conf.level) / 2),
                   "Upper" =  apply(mcmc, 1, quantile, probs = 1 - (1 - conf.level) / 2),
                   # Cumulative effect
                   "cEffect" = mean(apply(mcmc, 2, sum)),
                   "cLower" = quantile(apply(mcmc, 2, sum), (1 - conf.level) / 2),
                   "cUpper" = quantile(apply(mcmc, 2, sum), 1 - (1 - conf.level) / 2)))
      expDLM[[exposure]]$CW <- # Critical window where CI does not contain zero
        (expDLM[[exposure]]$Lower > 0 | expDLM[[exposure]]$Upper < 0)
    }
  }
  
  if (keep.mcmc) {
    return(expDLM)
  } else {
    return(do.call(rbind.data.frame, expDLM))
  }
}