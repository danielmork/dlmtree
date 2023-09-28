#' sim.tdlmm.zinb
#'
#' @description Simulation scenarios to accompany TDLM/TDLMM for count data
#'
#' @param sim integer (1-2) specifying simulation scenario
#' @param ctnum number of counties, 
#' @param week Weeks for each county, must be between 1-561
#' @param expList Named list of exposure data
#' @param data_zinb covariate data
#' @param zi_mean Mean of at-risk coefficients to control the zero-inflation of the count outcome
#' @param effect.size Parameter to control the magnitude of the count outcome
#' @param r Dispersion parameter
#'
#' @return A simulated dataset with true parameter values & zero-inflation information
#' @export
#'
sim.tdlmm.zinb <- function(sim = 1,
                           ctnum = 20,
                           week = 561,
                           expList = NULL,
                           data_zinb = NULL,
                           zi_mean = 0,
                           effect.size = 0.1,
                           r = 1) 
{
  if (!(sim %in% 1:2))
    stop("`sim` must be an integer from 1-2")

  if (week > length(which(data_zinb$fipscoor == unique(data_zinb$fipscoor)[1]))){
    stop(paste0("Week must be less than ", length(which(data_zinb$fipscoor == unique(data_zinb$fipscoor)[1]))))
  }

  # Exposure setup
  Lags <- ncol(expList[[1]])              # Lags: 40 weeks
  n <- week * ctnum
  n.samp <- min(nrow(expList[[1]]), n)    # minimum between {observation and the number of required sampling}

  # Sample from each county
  fips <- sample(unique(data_zinb$fipscoor), ctnum, replace = FALSE)
  
  idx = c()
  for(code in fips){
    fips_idx = sample(which(data_zinb$fipscoor == code), week, replace = FALSE)
    idx = c(idx, fips_idx)
  }
  
  # Data setup
  data_zinb <- data_zinb[idx, ]
  data_zinb$fipscoor <- droplevels(data_zinb$fipscoor)

  # Model matrix construction
  data_zi <- model.matrix(~ fipscoor - 1, data = data_zinb)
  data_nb <- model.matrix(~ fipscoor + month + YOC, data = data_zinb)

  # Sim 1: Zero-inflated negative binomial with single exposure DLM
  if (sim == 1) {
    # center/scale exposure data
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))

    # generate random starting time
    start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
    eff1 <- rep(0, Lags)                        # eff1 = a vector of 37 zeros
    eff1[start.time1:(start.time1 + 7)] <- 1    # Set 7 weeks of eff1 as 1

    # Start a y vector as 0.
    y <- rep(0, n.samp)

    # At-risk regression coefficients
    beta1 <- rnorm(ncol(data_zi), mean = zi_mean, sd = 1) # Sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)       # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          

    pi <- 1 / (1 + exp(-eta1))               # Inverse-logit
    w <- rbinom(n.samp, 1, pi)           # Zero-inflation indicator variable
    nStar <- n.samp - sum(w)                       # Proportion of structural zeros

    # NB regression coefficients
    beta2 <- rnorm(ncol(data_nb))                             # Sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 0, 1:ncol(data_nb)] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Update f to use only at-risk observations
    truth1Sums <- exposures[[1]][w == 0, ] %*% eff1           # A single exposure (n* x 37) %*% eff1 (37 x 1) = (n* x 1)
    f <- truth1Sums
    eta2_dlm <- effect.size * (eta2 + f)
    psi <- 1 / (1 + exp(-eta2_dlm))            # Probability of success in negative binomial

    # Draw y with eta1, eta2, and f
    y[w == 0] <- rnbinom(nStar, r, prob = 1 - psi)

    # Control the magnitude of the count outcome
    f <- f * effect.size
    margDLM1 <- eff1 * effect.size

    # Zero-inflation information
    zeroStr <- length(y[y == 0 & w == 1])/n   # y = 0 & classified as structural
    zeroAr <- length(y[y == 0 & w == 0])/n    # y = 0 & classified as at-risk
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

  # Sim 2: ZINB with 2 exposures with/without interaction
  if(sim == 2){
    # DLM Structure: Choose window of susceptibility for PM2.5 and temperature
    # Centering and Scaling
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))

    # Starting points & effects
    start.time1 <- sample(1:(Lags - 7), 1)
    start.time2 <- sample(1:(Lags - 7), 1)

    eff1 <- eff2 <- rep(0, Lags)
    eff1[start.time1:(start.time1 + 7)] <- 1
    eff2[start.time2:(start.time2 + 7)] <- 1

    # Start a y vector as 0.
    y <- rep(0, n.samp)

    # At-risk regression coefficients
    beta1 <- rnorm(ncol(data_zi), mean = zi_mean, sd = 1)   # Sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)         # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          

    pi <- 1 / (1 + exp(-eta1))           # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, pi)           # Zero-inflation indicator variable
    nStar <- n.samp - sum(w)                 

    # NB regression coefficients
    beta2 <- rnorm(ncol(data_nb))                             # Sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 0, 1:ncol(data_nb)] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Compute f with interaction
    truth1Sums <- exposures[[1]][w == 0, ] %*% eff1
    truth2Sums <- exposures[[2]][w == 0, ] %*% eff2
    int.effect = 0.025

    f <- truth1Sums + (int.effect * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect

    # NB with fixed and exposure effect
    eta2_dlm <- effect.size * (eta2 + f)
    psi <- 1 / (1 + exp(-eta2_dlm))            # Probability of success in negative binomial

    # Draw y with eta1, eta2, and f
    y[w == 0] <- rnbinom(nStar, r, prob = 1 - psi)

    # Calculate marginalized effects
    truthInt <- outer(eff1, eff2) * int.effect * effect.size
    margDLM1 <- eff1 * effect.size + rowSums(truthInt) * mean(exposures[[2]])
    margDLM2 <- colSums(truthInt) * mean(exposures[[1]])

    # Zero-inflation information
    zeroStr <- length(y[y == 0 & w == 1])/n   # y = 0 & classified as structural
    zeroAr <- length(y[y == 0 & w == 0])/n    # y = 0 & classified as at-risk
    nonzero <- length(y[y != 0])/n            # y is not zero hence structural
    zeroProp <- length(y[y == 0])/n           # Proportion of zeros

    return(list("data" = cbind.data.frame(y, data_zinb),
                "exposures" = exposures,
                "eff1" = eff1, "eff2" = eff2,
                "outer" = outer(eff1, eff2),
                "truthInt" = truthInt,
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
