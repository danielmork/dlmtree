#' sim.tdlmm
#'
#' @title Creates simulated data for TDLM & TDLMM
#' @description Method for creating simulated data for TDLM & TDLMM
#'
#' @param sim character (A - F) specifying simulation scenario
#' @param error positive scalar specifying error variance for Gaussian response
#' @param mean.p scalar between zero and one specifying mean probability for simulation scenario A
#' @param prop.active proportion of active exposures for simulation scenario C
#' @param n sample size for simulation
#' @param expList named list of exposure data
#' @param r dispersion parameter of negative binomial distribution
#'
#' @details Simulation scenarios:
#' - Scenario A: Binary response with single exposure effect
#' - Scenario B: Continuous response with main effect of PM2.5 and interaction
#' - Scenario C: Continuous response to test exposure selection using exposure
#' - Scenario D: Continuous response to test exposure selection using one exposure of main effect and two interaction effects among four exposures
#' - Scenario E: Zero-inflated count response with single exposure effect
#' - Scenario F: Zero-inflated count response with single exposure effect with main effect of PM2.5 and interaction
#' @md
#'
#' @examples
#' sim.tdlmm(sim = "A", mean.p = 0.5, n = 1000)
#'
#' @returns Simulated data and true parameters
#' @export
#'
sim.tdlmm <- function(sim = "A",
                      n = 5000,             # Scn A, B, C, D
                      error = 10,           # Scn B, C, D
                      mean.p = 0.5,         # Scn A
                      # n.exp = 25,         # Not in use (number of exposures for simulation scenarios three and four)
                      prop.active = 0.05,   # Scn C
                      expList = NULL,
                      r = 1)
{
  if (!(sim %in% LETTERS[1:6])) {
    stop("`sim` must be an character from A - F")
  }

  # TDLM/TDLMM
  if (sim %in% LETTERS[1:4]) {
    if (is.null(expList)) {
      data(exposureCov, envir = environment()) # Covariace matrix of exposure data
      cholexp     <- t(chol(exposureCov))
      exposureDat <- t(sapply(1:n, function(i) cholexp %*% rnorm(ncol(cholexp))))
      exposureDat <- exposureDat - min(exposureDat)
      expList     <- list("e1" = exposureDat[,1:37], 
                          "e2" = exposureDat[,38:74],
                          "e3" = exposureDat[,75:111],
                          "e4" = exposureDat[,112:148],
                          "e5" = exposureDat[,149:185])
    }

    Lags    <- ncol(expList[[1]])                                   # Lags: 37 weeks
    n.samp  <- min(nrow(expList[[1]]), n)                           # minimum between {observation and the number of required sampling}
    idx     <- sample(nrow(expList[[1]]), size = n.samp)            # sample n.samp# from the entire sample n. -> Used as indices
    data    <- cbind(matrix(rnorm(5 * n.samp), n.samp, 5),          # 5 guassian predictors (n.samp x 5)
                     matrix(rbinom(5 * n.samp, 1, .5), n.samp, 5))  # 5 binary predictors (n.samp x 5)       # => data: (n.samp x 10)
    
    colnames(data) <- c(paste0("c", 1:5), paste0("b", 1:5))         # Name the columns with c1 - c5 and b1 - b5
    params  <- rnorm(10)                                            # Sample true beta from a standard normal
    c       <- c(data[,1:10] %*% params)                            # Compute c: zT * gamma (nx10) x (10x1) = (nx1)


    # Sim A: binary response with single exposure DLM
    if (sim == "A") {
      exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,])) # center/scale exposure data

      start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
      eff1        <- rep(0, Lags)                 # eff1 = a vector of 37 zeros
      eff1[start.time1:(start.time1 + 7)] <- 1    # Set 7 weeks of eff1 as 1
      truth1Sums  <- exposures[[1]] %*% eff1     
      f           <- truth1Sums

      # generate output and set alpha for mean.p
      pp  <- (-0.1) * (c + f) 
      p   <- 1 / (1 + exp(pp - mean(pp) + log(mean.p^(-1) - 1)))  # Calculate a vector of Bernoulli probability, p
      f   <- 0.1 * f # return the true vector of f

      # draw y
      y   <- rbinom(n.samp, 1, p)
      margDLM1 <- eff1 * 0.1
      return(list("dat" = cbind.data.frame(y, data), 
                  "params" = params * 0.1, 
                  "exposures" = exposures, 
                  "start.time1" = start.time1,
                  "margDLM" = margDLM1, 
                  "c" = c, 
                  "f" = f, 
                  "p" = p))
    }


    # Sim B: continuous response with main effect of PM2.5 and interaction between PM2.5 and NO2
    if (sim == "B") {
      exposures     <- lapply(expList, function(i) i[idx,] / IQR(i[idx,]))
      start.time1   <- sample(1:(Lags - 7), 1)
      start.time2   <- sample(1:(Lags - 7), 1)

      eff1 <- eff2 <- rep(0, Lags)
      eff1[start.time1:(start.time1 + 7)] <- 1
      eff2[start.time2:(start.time2 + 7)] <- 1
      truth1Sums <- exposures[[1]] %*% eff1 
      truth2Sums <- exposures[[2]] %*% eff2 

      f <- truth1Sums + (0.025 * truth1Sums * truth2Sums) # The critical window + overlap interaction
      effect.size <- 1 / sd(f)
      f <- f * effect.size

      # draw y
      y         <- c + f + rnorm(n.samp, sd = sqrt(error))
      truthInt  <- outer(eff1, eff2) * 0.025 * effect.size 
      margDLM1  <- eff1 * effect.size + rowSums(truthInt) * mean(exposures[[2]])
      margDLM2  <- colSums(truthInt) * mean(exposures[[1]])

      return(list("dat" = cbind.data.frame(y, data),
                  "exposures" = exposures, 
                  "params" = params,
                  "start.time1" = start.time1, 
                  "start.time2" = start.time2,
                  "margDLM1" = margDLM1, 
                  "margDLM2" = margDLM2,
                  "effectsize" = effect.size,
                  "c" = c, 
                  "f" = f))
    }

    # Create exposures for scenarios C & D
    S11       <- exp(-toeplitz(0:(36)) * 0.7)
    a         <- 0.5
    exp       <- list()
    exp[[1]]  <- t(sapply(1:n, function(i) t(chol(S11)) %*% rnorm(37)))
    for (i in 2:p) {
      exp[[i]] <- t(sapply(1:n, function(j) {
        t(chol(S11)) %*% rnorm(37, sd = sqrt(1 - a^2)) + a * exp[[i - 1]][j,] }))
    }
    names(exp) <- paste0("e", 1:p)

    # Sim C: continuous response to test exposure selection using exposure
    # main effects
    if (sim == "C") {
      # Draw active exposures and assign DLM effects
      active      <- sample.int(p, round(prop.active * p))
      active.dlm  <- list()
      f <- rep(0, n)
      for (i in active) {
        start           <- sample.int(30, 1)
        active.dlm[[i]] <- start:(start + 7)
        f               <- f + rowSums(exp[[i]][,start:(start + 7)])
      }
      f <- scale(f)

      # Create response
      y <- c + f + rnorm(n, sd = sqrt(error))

      return(list("dat" =  cbind.data.frame(y, data), 
                  "params" = params,
                  "exposures" = exp,
                  "active" = active,
                  "active.dlm" = active.dlm,
                  "f" = f,
                  "c" = c))
    }

    # Sim 4: continuous response to test exposure selection using one exposure
    # main effect and two interaction effects
    if (sim == "D") {
      # Draw active exposures and assign DLM effects
      active      <- sort(sample.int(p, 5))
      active.int  <- c(paste0(names(exp)[active[2]], "-", names(exp)[active[3]]),
                       paste0(names(exp)[active[4]], "-", names(exp)[active[5]]))

      start1 <- sample.int(30, 1)
      start2 <- sample.int(30, 1)
      start3 <- sample.int(30, 1)
      start4 <- sample.int(30, 1)
      start5 <- sample.int(30, 1)
      f <- rowSums(exp[[active[1]]][,start1:(start1 + 7)]) +
        0.2 * rowSums(exp[[active[2]]][,start2:(start2 + 7)]) *
        rowSums(exp[[active[3]]][,start3:(start3 + 7)]) +
        0.2 * rowSums(exp[[active[4]]][,start4:(start4 + 7)]) *
        rowSums(exp[[active[5]]][,start5:(start5 + 7)])
      f <- scale(f)

      # Create response
      y <- c + f + rnorm(n, sd = sqrt(error))

      return(list("dat" =  cbind.data.frame(y, data), 
                  "params" = params,
                  "exposures" = exp,
                  "active" = active,
                  "active.int" = active.int,
                  "f" = f, 
                  "c" = c))
    }
  } else { 
    # Read in data containing exposure & covariates
    data(zinbCo, envir = environment())
    week    <- 52
    ctnum   <- 64

    # Time-series covariates
    birthCov <- zinbCo[, c("fipscoor", "month", "YOC")]

    # Exposure data
    if (is.null(expList)) {
      expList <- list(
        e1 = as.matrix(zinbCo[, grep("^cmaq_pm25_", names(zinbCo), value = TRUE)]),
        e2 = as.matrix(zinbCo[, grep("^tmmx_", names(zinbCo), value = TRUE)])
      )
    } else {
      if(length(expList) < 2){
        stop("The length of an input list 'expList' must be equal to or more than two.")
      }
    }

    # Setup
    Lags    <- ncol(expList[[1]])           # lags: 40 weeks
    n       <- week * ctnum                 # sample size
    n.samp  <- min(nrow(expList[[1]]), n)   # minimum between {observation and the number of required sampling}

    # Sample counties
    fips <- sample(unique(birthCov$fipscoor), ctnum, replace = FALSE)
    
    # Sample weeks from each county
    idx <- c()
    for (code in fips) {
      fips.idx  <- sample(which(birthCov$fipscoor == code), week, replace = FALSE)
      idx       <- c(idx, fips.idx)
    }
    
    # Data setup
    birthCov           <- birthCov[idx, ]
    birthCov$fipscoor  <- droplevels(birthCov$fipscoor)

    # Model matrix construction for ZI & NB model
    data.zi <- model.matrix(~ fipscoor - 1, data = birthCov)
    data.nb <- model.matrix(~ fipscoor + month + YOC, data = birthCov)

    # Effect size
    effect.size <- 0.1

    # Sim E: A single exposure DLM
    if (sim == "E") {
      exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,])) # center/scale exposure data

      # Sample starting time lag
      start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
      eff1        <- rep(0, Lags)                 # eff1 = a vector of 37 zeros
      eff1[start.time1:(start.time1 + 7)] <- 1    # set 8 weeks of eff1 as 1

      # Start y vector as 0
      y <- rep(0, n.samp)

      # ZI regression coefficients
      beta1 <- rnorm(ncol(data.zi), mean = 0, sd = 1)   # sample true beta1 from a standard normal
      eta1  <- c(data.zi[, 1:ncol(data.zi)] %*% beta1)  # ZI model fixed effect    

      pi    <- 1 / (1 + exp(-eta1))                     # inverse-logit
      w     <- rbinom(n.samp, 1, pi)                    # ZI indicator variable
      nStar <- n.samp - sum(w)                          # number of NB observations

      # NB regression coefficients
      beta2 <- rnorm(ncol(data.nb))                            # sample true beta2 from a standard normal
      eta2  <- c(data.nb[w == 0, 1:ncol(data.nb)] %*% beta2)   # NB model fixed effect

      # Update f to use only NB observations
      f         <- exposures[[1]][w == 0, ] %*% eff1    # NB model exposure effect
      eta2.dlm  <- effect.size * (eta2 + f)             # scaled NB model exposure effect
      psi       <- 1 / (1 + exp(-eta2.dlm))             # probability of success in NB model

      # Draw y with eta1, eta2, and f
      y[w == 0] <- rnbinom(nStar, r, prob = 1 - psi) 

      # Scale the magnitude of the count outcome for returning the outcome
      f         <- f * effect.size
      margDLM1  <- eff1 * effect.size

      # ZI information
      zeroStr   <- length(y[y == 0 & w == 1])/n   # y = 0 & classified as ZI
      zeroAr    <- length(y[y == 0 & w == 0])/n   # y = 0 & classified as NB
      nonzero   <- length(y[y != 0])/n            # y is not zero
      zeroProp  <- length(y[y == 0])/n            # Proportion of zeros

      # Return the result
      return(list("data" = cbind.data.frame(y, birthCov),
                  "exposures" = exposures, 
                  "start.time1" = start.time1,
                  "eff1" = eff1,
                  "margDLM1" = margDLM1, 
                  "eta1" = eta1, 
                  "eta2" = eta2,
                  "b1" = beta1, 
                  "b2" = beta2 * effect.size,
                  "pi" = pi, 
                  "psi" = psi, 
                  "r" = r, 
                  "f" = f, 
                  "w" = w,
                  "zeroStr" = zeroStr,
                  "zeroAr" = zeroAr,
                  "nonzero" = nonzero,
                  "zeroProportion" = zeroProp))
    }

    # Sim F: Two exposures with interaction
    if (sim == "F") {
      exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))  # centering and scaling

      start.time1 <- sample(1:(Lags - 7), 1)    # e1 starting time lag
      start.time2 <- sample(1:(Lags - 7), 1)    # e2 starting time lag

      eff1 <- eff2 <- rep(0, Lags)              
      eff1[start.time1:(start.time1 + 7)] <- 1  # e1 lag effect
      eff2[start.time2:(start.time2 + 7)] <- 1  # e2 lag effect

      # Start y vector as 0.
      y <- rep(0, n.samp)

      # ZI regression coefficients
      beta1 <- rnorm(ncol(data.zi), mean = 0, sd = 1)    # Sample true beta1 from a standard normal
      eta1  <- c(data.zi[, 1:ncol(data.zi)] %*% beta1)   # ZI model fixed effect

      pi    <- 1 / (1 + exp(-eta1))    # 1 - P(ZI zero)
      w     <- rbinom(n.samp, 1, pi)   # ZI indicator variable
      nStar <- n.samp - sum(w)                 

      # NB regression coefficients
      beta2 <- rnorm(ncol(data.nb))                             # Sample true beta2 from a standard normal
      eta2  <- c(data.nb[w == 0, 1:ncol(data.nb)] %*% beta2)    # NB model fixed effect

      # Compute f with interaction
      truth1Sums <- exposures[[1]][w == 0, ] %*% eff1     # e1 exposure effect
      truth2Sums <- exposures[[2]][w == 0, ] %*% eff2     # e2 exposure effect
      int.effect <- 0.025

      f <- truth1Sums + (int.effect * truth1Sums * truth2Sums)  # e1 main effect + e1xe2 interaction effect

      # NB with fixed and exposure effect
      eta2.dlm  <- effect.size * (eta2 + f)
      psi       <- 1 / (1 + exp(-eta2.dlm))  # Probability of success in negative binomial

      # Draw y with eta1, eta2, and f
      y[w == 0] <- rnbinom(nStar, r, prob = 1 - psi)

      # Calculate marginalized effects
      truthInt <- outer(eff1, eff2) * int.effect * effect.size
      margDLM1 <- eff1 * effect.size + rowSums(truthInt) * mean(exposures[[2]])
      margDLM2 <- colSums(truthInt) * mean(exposures[[1]])

      # ZI information
      zeroStr   <- length(y[y == 0 & w == 1])/n  # y = 0 & classified as ZI
      zeroAr    <- length(y[y == 0 & w == 0])/n  # y = 0 & classified as NB
      nonzero   <- length(y[y != 0])/n           # y is not zero
      zeroProp  <- length(y[y == 0])/n           # Proportion of zeros

      return(list("data" = cbind.data.frame(y, birthCov),
                  "exposures" = exposures,
                  "eff1" = eff1, 
                  "eff2" = eff2,
                  "start.time1" = start.time1, 
                  "start.time2" = start.time2,
                  "margDLM1" = margDLM1,
                  "margDLM2" = margDLM2,
                  "f" = f, 
                  "w" = w, 
                  "psi" = psi, 
                  "r" = r,
                  "eta1" = eta1, 
                  "eta2" = eta2.dlm,
                  "b1" = beta1, 
                  "b2" = beta2 * effect.size,
                  "zeroStr" = zeroStr,
                  "zeroAr" = zeroAr,
                  "nonzero" = nonzero,
                  "zeroProportion" = zeroProp))
      }
    }
}
