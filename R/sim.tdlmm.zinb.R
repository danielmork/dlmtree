#' sim.tdlmm.zinb
#'
#' @description Simulation scenarios to accompany TDLM/TDLMM for count data
#'
#' @param sim integer (1-2) specifying simulation scenario
#' @param ctnum number of counties
#' @param week number of weeks for each county, must be between 1-561
#' @param expList named list of exposure data
#' @param data_zinb covariate data
#' @param effect.size the effect size of the window of susceptibility, resulting in controlling the magnitude of the count outcome
#' @param r dispersion parameter
#'
#' @return A simulated dataset, true parameter values, zero/non-zero proportions
#' @export
#'
sim.tdlmm.zinb <- function(sim = 1,
                           ctnum = 20,
                           week = 561,
                           expList = NULL,
                           data_zinb = NULL,
                           effect.size = 0.1,
                           r = 1) 
{
  if (!(sim %in% 1:2))
    stop("`sim` must be an integer from 1-2")

  if (week > length(which(data_zinb$fipscoor == unique(data_zinb$fipscoor)[1]))){
    stop(paste0("Week must be less than ", length(which(data_zinb$fipscoor == unique(data_zinb$fipscoor)[1]))))
  }

  # Setup
  Lags <- ncol(expList[[1]])              # lags: 40 weeks
  n <- week * ctnum                       # sample size
  n.samp <- min(nrow(expList[[1]]), n)    # minimum between {observation and the number of required sampling}

  # Sample counties
  fips <- sample(unique(data_zinb$fipscoor), ctnum, replace = FALSE)
  
  # Sample weeks from each county
  idx = c()
  for(code in fips){
    fips_idx = sample(which(data_zinb$fipscoor == code), week, replace = FALSE)
    idx = c(idx, fips_idx)
  }
  
  # Data setup
  data_zinb <- data_zinb[idx, ]
  data_zinb$fipscoor <- droplevels(data_zinb$fipscoor)

  # Model matrix construction for ZI & NB model
  data_zi <- model.matrix(~ fipscoor - 1, data = data_zinb)
  data_nb <- model.matrix(~ fipscoor + month + YOC, data = data_zinb)

  # Sim 1: A single exposure DLM
  if (sim == 1) {
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,])) # center/scale exposure data

    # Sample starting time lag
    start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
    eff1 <- rep(0, Lags)                        # eff1 = a vector of 37 zeros
    eff1[start.time1:(start.time1 + 7)] <- 1    # set 8 weeks of eff1 as 1

    # Start y vector as 0
    y <- rep(0, n.samp)

    # ZI regression coefficients
    beta1 <- rnorm(ncol(data_zi), mean = 0, sd = 1)   # sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)   # ZI model fixed effect    

    pi <- 1 / (1 + exp(-eta1))                        # inverse-logit
    w <- rbinom(n.samp, 1, pi)                        # ZI indicator variable
    nStar <- n.samp - sum(w)                          # number of NB observations

    # NB regression coefficients
    beta2 <- rnorm(ncol(data_nb))                             # sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 0, 1:ncol(data_nb)] %*% beta2)     # NB model fixed effect

    # Update f to use only NB observations
    f <- exposures[[1]][w == 0, ] %*% eff1     # NB model exposure effect
    eta2_dlm <- effect.size * (eta2 + f)       # scaled NB model exposure effect
    psi <- 1 / (1 + exp(-eta2_dlm))            # probability of success in NB model

    # Draw y with eta1, eta2, and f
    y[w == 0] <- rnbinom(nStar, r, prob = 1 - psi) 

    # Scale the magnitude of the count outcome for returning the outcome
    f <- f * effect.size
    margDLM1 <- eff1 * effect.size

    # ZI information
    zeroStr <- length(y[y == 0 & w == 1])/n   # y = 0 & classified as ZI
    zeroAr <- length(y[y == 0 & w == 0])/n    # y = 0 & classified as NB
    nonzero <- length(y[y != 0])/n            # y is not zero
    zeroProp <- length(y[y == 0])/n           # Proportion of zeros

    # Return the result
    return(list("data" = cbind.data.frame(y, data_zinb),
                "exposures" = exposures, 
                "start.time1" = start.time1,
                "eff1" = eff1,
                "margDLM1" = margDLM1, 
                "eta1" = eta1, "eta2" = eta2,
                "b1" = beta1, "b2" = beta2 * effect.size,
                "pi" = pi, "psi" = psi, "r" = r, "f" = f, "w" = w,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp))
  }

  # Sim 2: Two exposures with interaction
  if(sim == 2){
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))     # centering and scaling

    start.time1 <- sample(1:(Lags - 7), 1)    # e1 starting time lag
    start.time2 <- sample(1:(Lags - 7), 1)    # e2 starting time lag

    eff1 <- eff2 <- rep(0, Lags)              
    eff1[start.time1:(start.time1 + 7)] <- 1  # e1 lag effect
    eff2[start.time2:(start.time2 + 7)] <- 1  # e2 lag effect

    # Start y vector as 0.
    y <- rep(0, n.samp)

    # ZI regression coefficients
    beta1 <- rnorm(ncol(data_zi), mean = 0, sd = 1)         # Sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)         # ZI model fixed effect

    pi <- 1 / (1 + exp(-eta1))           # 1 - P(ZI zero)
    w <- rbinom(n.samp, 1, pi)           # ZI indicator variable
    nStar <- n.samp - sum(w)                 

    # NB regression coefficients
    beta2 <- rnorm(ncol(data_nb))                             # Sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 0, 1:ncol(data_nb)] %*% beta2)     # NB model fixed effect

    # Compute f with interaction
    truth1Sums <- exposures[[1]][w == 0, ] %*% eff1     # e1 exposure effect
    truth2Sums <- exposures[[2]][w == 0, ] %*% eff2     # e2 exposure effect
    int.effect = 0.025

    f <- truth1Sums + (int.effect * truth1Sums * truth2Sums)  # e1 main effect + e1xe2 interaction effect

    # NB with fixed and exposure effect
    eta2_dlm <- effect.size * (eta2 + f)
    psi <- 1 / (1 + exp(-eta2_dlm))            # Probability of success in negative binomial

    # Draw y with eta1, eta2, and f
    y[w == 0] <- rnbinom(nStar, r, prob = 1 - psi)

    # Calculate marginalized effects
    truthInt <- outer(eff1, eff2) * int.effect * effect.size
    margDLM1 <- eff1 * effect.size + rowSums(truthInt) * mean(exposures[[2]])
    margDLM2 <- colSums(truthInt) * mean(exposures[[1]])

    # ZI information
    zeroStr <- length(y[y == 0 & w == 1])/n   # y = 0 & classified as ZI
    zeroAr <- length(y[y == 0 & w == 0])/n    # y = 0 & classified as NB
    nonzero <- length(y[y != 0])/n            # y is not zero
    zeroProp <- length(y[y == 0])/n           # Proportion of zeros

    return(list("data" = cbind.data.frame(y, data_zinb),
                "exposures" = exposures,
                "eff1" = eff1, "eff2" = eff2,
                "start.time1" = start.time1, "start.time2" = start.time2,
                "margDLM1" = margDLM1, "margDLM2" = margDLM2,
                "f" = f, "w" = w, "psi" = psi, "r" = r,
                "eta1" = eta1, "eta2" = eta2_dlm,
                "b1" = beta1, "b2" = beta2 * effect.size,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp))
  }
}
