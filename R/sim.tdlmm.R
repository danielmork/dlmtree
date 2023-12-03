#' sim.tdlmm
#'
#' @description Simulation scenarios to accompany TDLM/TDLMM
#'
#' @param sim integer (1-6) specifying simulation scenario
#' @param error positive scalar specifying error variance
#' @param mean.p scalar between zero and one specifying mean probability
#' for simulation scenario one
#' @param n.exp number of exposures for simulation scenarios three and four
#' @param prop.active proportion of active exposures for simulation scenario
#' three
#' @param n sample size for simulation
#' @param expList named list of exposure data
#' @param ctnum number of counties
#' @param week number of weeks for each county, must be between 1-561
#' @param expList named list of exposure data
#' @param data_zinb covariate data
#' @param effect.size the effect size of the window of susceptibility, resulting in controlling the magnitude of the count outcome
#' @param r dispersion parameter
#'
#' @return Simulated data and true parameters
#' @export
#'
sim.tdlmm <- function(sim = 1,
                      error = 10,           # Scn 2, 3, 4
                      mean.p = 0.5,         # Scn 1
                      n.exp = 25,           # Not in use
                      prop.active = 0.05,   # Scn 3
                      n = 5000,             # Scn 1, 2, 3, 4
                      expList = NULL,
                      # ZINB
                      ctnum = 20,           # n = ctnum * week for zinb (time series)
                      week = 561,
                      data_zinb = NULL,
                      r = 1)
{
  if (!(sim %in% 1:6))
    stop("`sim` must be an integer from 1-4")

  # TDLM/TDLMM (Dan)
  if(sim %in% 1:4){
    if (is.null(expList)) {
      data(exposureCov) # Covariace matrix of exposure data
      cholexp <- t(chol(exposureCov))
      exposureDat <- t(sapply(1:n, function(i) cholexp %*% rnorm(ncol(cholexp))))
      exposureDat <- exposureDat - min(exposureDat)
      expList <- list("e1" = exposureDat[,1:37], # Each exposure has 37 weeks
                      "e2" = exposureDat[,38:74],
                      "e3" = exposureDat[,75:111],
                      "e4" = exposureDat[,112:148],
                      "e5" = exposureDat[,149:185])
    }

    Lags <- ncol(expList[[1]])              # Lags: 37 weeks
    n.samp <- min(nrow(expList[[1]]), n)    # minimum between {observation and the number of required sampling}
    idx <- sample(nrow(expList[[1]]), size = n.samp) # sample n.samp# from the entire sample n. -> Used as indices
    data <- cbind(matrix(rnorm(5 * n.samp), n.samp, 5),           # 5 guassian predictors (n.samp x 5)
                  matrix(rbinom(5 * n.samp, 1, .5), n.samp, 5))   # 5 binary predictors (n.samp x 5)       # => data: (n.samp x 10)
    colnames(data) <- c(paste0("c", 1:5), paste0("b", 1:5))       # Name the columns with c1 - c5 and b1 - b5
                                                                  # c for continuous, b for binary
    params <- rnorm(10)                                           # Sample true beta from a standard normal
    c <- c(data[,1:10] %*% params)                                # Compute c: zT * gamma (nx10) x (10x1) = (nx1)


    # Sim 1: binary response with single exposure DLM
    if (sim == 1) {
      # center/scale exposure data
      exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))
      # generate random starting time
      start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
      eff1 <- rep(0, Lags)                        # eff1 = a vector of 37 zeros
      eff1[start.time1:(start.time1 + 7)] <- 1    # Set 7 weeks of eff1 as 1
      truth1Sums <- exposures[[1]] %*% eff1       # A single exposure (n x 37) %*% eff1 (37 x 1) = (n x 1)
      f <- truth1Sums # (n x 1)
      # generate output and set alpha for mean.p
      pp <- (-0.1) * (c + f) # (n x 1)
      p <- 1 / (1 + exp(pp - mean(pp) + log(mean.p^(-1) - 1)))  # Calculate a vector of Bernoulli probability, p
      f <- 0.1 * f # return the true vector of f
      # draw y
      y <- rbinom(n.samp, 1, p)
      margDLM1 <- eff1 * 0.1
      return(list("dat" = cbind.data.frame(y, data), 
                  "params" = params * 0.1, # return the true coefficient of 10 parameters
                  "exposures" = exposures, 
                  "start.time1" = start.time1,
                  "margDLM" = margDLM1, 
                  "c" = c, "f" = f, "p" = p))
    }


    # Sim 2: continuous response with main effect of PM2.5 and interaction
    # between PM2.5 and NO2
    if (sim == 2) {
      exposures <- lapply(expList, function(i) i[idx,] / IQR(i[idx,]))
      start.time1 <- sample(1:(Lags - 7), 1)
      start.time2 <- sample(1:(Lags - 7), 1)
      eff1 <- eff2 <- rep(0, Lags)
      eff1[start.time1:(start.time1 + 7)] <- 1
      eff2[start.time2:(start.time2 + 7)] <- 1
      truth1Sums <- exposures[[1]] %*% eff1 # (nxT)(Tx1) = (nx1) # Sum of effects on the critical window
      truth2Sums <- exposures[[2]] %*% eff2 # (nxT)(Tx1) = (nx1) # Sum of effects on the critical window
      f <- truth1Sums + (0.025 * truth1Sums * truth2Sums) # The critical window + overlap interaction
      effect.size <- 1 / sd(f)
      f <- f * effect.size
      # draw y
      y <- c + f + rnorm(n.samp, sd = sqrt(error))
      truthInt <- outer(eff1, eff2) * 0.025 * effect.size #True interaction
      margDLM1 <- eff1 * effect.size + rowSums(truthInt) * mean(exposures[[2]])# effect size * (eff1 + 0.025*E[e2])
      margDLM2 <- colSums(truthInt) * mean(exposures[[1]])# effect size * E[e1]
      return(list("dat" = cbind.data.frame(y, data),
                  "exposures" = exposures, "params" = params,
                  "start.time1" = start.time1, "start.time2" = start.time2,
                  "margDLM1" = margDLM1, "margDLM2" = margDLM2,
                  "effectsize" = effect.size,
                  "c" = c, "f" = f))
    }

    # create exposures for scenarios three and four
    S11 <- exp(-toeplitz(0:(36)) * 0.7)
    a <- 0.5
    exp <- list()
    exp[[1]] <- t(sapply(1:n, function(i) t(chol(S11)) %*% rnorm(37)))
    for (i in 2:p) {
      exp[[i]] <- t(sapply(1:n, function(j) {
        t(chol(S11)) %*% rnorm(37, sd = sqrt(1 - a^2)) + a * exp[[i - 1]][j,] }))
    }
    names(exp) <- paste0("e", 1:p)


    # Sim 3: continuous response to test exposure selection using exposure
    # main effects
    if (sim == 3) {
      # Draw active exposures and assign DLM effects
      active <- sample.int(p, round(prop.active * p))
      active.dlm <- list()
      f <- rep(0, n)
      for (i in active) {
        start <- sample.int(30, 1)
        active.dlm[[i]] <- start:(start + 7)
        f <- f + rowSums(exp[[i]][,start:(start + 7)])
      }
      f <- scale(f)

      # Create response
      y <- c + f + rnorm(n, sd = sqrt(error))

      return(list("dat" =  cbind.data.frame(y, data), "params" = params,
                  "exposures" = exp,
                  "active" = active,
                  "active.dlm" = active.dlm,
                  "f" = f, "c" = c))
    }

    # Sim 4: continuous response to test exposure selection using one exposure
    # main effect and two interaction effects
    if (sim == 4) {
      # Draw active exposures and assign DLM effects
      active <- sort(sample.int(p, 5))
      active.int <- c(paste0(names(exp)[active[2]], "-", names(exp)[active[3]]),
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

      return(list("dat" =  cbind.data.frame(y, data), "params" = params,
                  "exposures" = exp,
                  "active" = active,
                  "active.int" = active.int,
                  "f" = f, "c" = c))
    }
  } else { 
    # ZINB (Seongwon)
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

    # Effect size
    effect.size = 0.1

    # Sim 5: A single exposure DLM
    if (sim == 5) {
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

    # Sim 6: Two exposures with interaction
    if(sim == 6){
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
}
