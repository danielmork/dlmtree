#' sim.tdlmm.zinb
#'
#' @description Simulation scenarios to accompany TDLM/TDLMM for count data
#'
#' @param sim integer (1-3) specifying simulation scenario
#' @param ctnum number of counties (or spatial units)
#' @param week Weeks for each county
#' @param expList Named list of exposure data
#' @param data_zinb data
#'
#' @return
#' @export
#'
sim.tdlmm.zinb <- function(sim = 1,
                            ctnum = 20, # Each county has 561 weeks
                            week = 561,
                            expList = NULL,
                            data_zinb = NULL,
                            rho = 0.95,
                            spTau = 1)
{
  if (!(sim %in% 1:3))
    stop("`sim` must be an integer from 1-3")

  if (week > length(which(data_zinb$fipscoor == unique(data_zinb$fipscoor)[1]))){
    stop(paste0("Week must be less than ", length(which(data_zinb$fipscoor == unique(data_zinb$fipscoor)[1]))))
  }

  # Exposure setup
  Lags <- ncol(expList[[1]])              # Lags: 40 weeks
  n <- week * ctnum
  n.samp <- min(nrow(expList[[1]]), n)    # minimum between {observation and the number of required sampling}

  if(sim == 3){
    # Collecting counties that are adjacent to each other
    start_cty <- sample(unique(data_zinb$fipscoor), 1)
    counties <- c(start_cty)
    i = 1

    while(length(counties) <= ctnum){
      start = start_cty
      new_cty <- suppressMessages(getAdjacentCounties(start))
      
      if(new_cty %in% counties){
        i = i + 1
              
        # In case sampling gets stuck in a county with one adjacency
        if(i %% 50 == 0){
          start_cty <- sample(unique(data_zinb$fipscoor), 1)
          counties = c(start_cty)
        } else {
          next
        }
      } else {
        counties <- c(counties, new_cty)
        start_cty <- new_cty
      }
    }

    fips <- counties[-1]
    W <- getAdjacencyMatrix(counties)
    D <- diag(rowSums(W))

  } else {
    fips <- sample(unique(data_zinb$fipscoor), ctnum, replace = FALSE)
  }

  # Sample from each county
  idx = c()
  for(code in fips){
    fips_idx = sample(which(data_zinb$fipscoor == code), week, replace = FALSE)
    idx = c(idx, fips_idx)
  }
  
  # Data setup
  data_zinb <- data_zinb[idx, ]
  data_zi <- model.matrix(~ fipscoor - 1, data = data_zinb)
  data_nb <- model.matrix(~ fipscoor + month + YOC, data = data_zinb)



  # Sim 1: Zero-inflated negative binomial with single exposure DLM (temporarily here as "p" is not defined)
  if (sim == 1) {
    # center/scale exposure data
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))
    # generate random starting time
    start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
    eff1 <- rep(0, Lags)                        # eff1 = a vector of 37 zeros
    eff1[start.time1:(start.time1 + 7)] <- 1    # Set 7 weeks of eff1 as 1

    # Start a y vector as 0.
    y <- rep(0, n.samp)

    # Need to redefine c (renamed to eta1 & eta2) as we have beta1 for binary component and beta2 for negbin component
    beta1 <- rnorm(ncol(data_zi)) # Sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)     # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          

    phi <- 1 / (1 + exp(-eta1))           # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, phi)           # "At-risk" indicator variable
    nStar <- sum(w)                       # Proportion of structural zeros

    # Sample beta2 for negbin component and use the variable w to calculate eta2 (At-risk observations only)
    beta2 <- rnorm(ncol(data_nb))  # Sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 1, 1:ncol(data_nb)] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Update f to use only at-risk
    truth1Sums <- exposures[[1]][w == 1, ] %*% eff1       # A single exposure (n* x 37) %*% eff1 (37 x 1) = (n* x 1)
    f <- truth1Sums # (n x 1)
    effect.size = 1 # 1/sd(f)
    f <- f * effect.size

    scale.size = 0.1
    eta2_dlm <- scale.size*(eta2 + f)
    pi <- 1 / (1 + exp(-eta2_dlm))            # Probability of success in negative binomial

    # draw y with eta1, eta2, and f
    r <- 1                                        # Dispersion parameter
    mu <- r * pi / (1 - pi)
    y[w == 1] <- rnbinom(nStar, r, mu = mu) 

    # return the true vector of f
    f <- scale.size * f
    margDLM1 <- eff1 * scale.size * effect.size

    zeroStr <- length(y[y == 0 & w == 0])/n # y = 0 & classified as structural
    zeroAr <- length(y[y == 0 & w != 0])/n  # y = 0 & classified as at-risk
    nonzero <- length(y[y != 0])/n # y is not zero hence structural
    zeroProp <- length(y[y == 0])/n  # Proportion of zeros: Note that this is heavily dependent on beta1 and beta2


    # return the result
    return(list("data" = cbind.data.frame(y, data_zinb),
                "exposures" = exposures, 
                "start.time1" = start.time1,
                "margDLM1" = margDLM1, 
                "eta1" = eta1, "eta2" = eta2,
                "b1" = beta1, "b2" = beta2 * scale.size,
                "phi" = phi, # At-risk probability
                "mu" = mu, "pi" = pi, "r" = r, "f" = f, "w" = w,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp))
  }

  # Sim 2: ZINB with 5 exposures with/without interaction
  if(sim == 2){
    # DLM Structure: Choose critical windows for PM2.5 and Temperature
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))
    start.time1 <- sample(1:(Lags - 7), 1)
    start.time2 <- sample(1:(Lags - 7), 1)
    eff1 <- eff2 <- rep(0, Lags)
    eff1[start.time1:(start.time1 + 7)] <- 1
    eff2[start.time2:(start.time2 + 7)] <- 1

    # Data generating process
    # Start a y vector as 0.
    y <- rep(0, n.samp)

    # Need to redefine c (renamed to eta1 & eta2) as we have beta1 for binary component and beta2 for negbin component
    beta1 <- rnorm(ncol(data_zi)) # Sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)     # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          

    phi <- 1 / (1 + exp(-eta1))           # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, phi)           # "At-risk" indicator variable
    nStar <- sum(w)                       # Proportion of at-risk zeros

    # Sample beta2 for negbin component and use the variable w to calculate eta2 (At-risk observations only)
    beta2 <- rnorm(ncol(data_nb))  # Sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 1, 1:ncol(data_nb)] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Compute f
    truth1Sums <- exposures[[1]][w == 1, ] %*% eff1
    truth2Sums <- exposures[[2]][w == 1, ] %*% eff2
    #f <- truth1Sums + (0.025 * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    int.effect = 0.025
    f <- truth1Sums + (int.effect * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    effect.size <- 1
    f <- f * effect.size

    scale.size = 0.1
    eta2_dlm <- scale.size * (eta2 + f)
    pi <- 1 / (1 + exp(-eta2_dlm))            # Probability of success in negative binomial

    # draw y with eta1, eta2, and f
    r <- 1                                  # Dispersion parameter
    mu <- r * pi / (1 - pi)
    y[w == 1] <- rnbinom(nStar, r, mu = mu)

    # Compute marginal effects
    # Scalar such that the variance of f is 1
    # f <- f * effect.size # Scaling 

    truthInt <- outer(eff1, eff2) * int.effect * effect.size * scale.size
    margDLM1 <- eff1 * effect.size * scale.size + rowSums(truthInt) * mean(exposures[[2]])
    margDLM2 <- colSums(truthInt) * mean(exposures[[1]])

    zeroStr <- length(y[y == 0 & w == 0])/n # y = 0 & classified as structural
    zeroAr <- length(y[y == 0 & w == 1])/n  # y = 0 & classified as at-risk
    nonzero <- length(y[y != 0])/n # y is not zero hence structural
    zeroProp <- length(y[y == 0])/n

    return(list("data" = cbind.data.frame(y, data_zinb),
                "exposures" = exposures,
                "eff1" = eff1,
                "eff2" = eff2,
                "outer" = outer(eff1, eff2),
                "truthInt" = truthInt,
                "start.time1" = start.time1, "start.time2" = start.time2,
                "margDLM1" = margDLM1, "margDLM2" = margDLM2,
                "c" = c, "f" = f, "w" = w, "mu" = mu, "phi" = phi, 
                "eta1" = eta1, "r" = r,
                "b1" = beta1,
                "b2" = beta2 * scale.size,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp))
  }

  # Sim 3: ZINB + TDLMM + Spatial random effect
  if(sim == 3){
    # CAR var-cov matrix
    Q <- D - rho*W
    Qinv <- spTau * chol2inv(chol(Q))

    # Compute the spatial edge nodes for faster computation
    spNodes = spatial_pairwise(W)

    # Compute variance-covariance matrix for phi (CAR model)
    sp_phi = c(mvrnorm(n = 1, rep(0, ctnum), Qinv)) # Simulate spatial effect

    # model matrix of spatial phi
    A = matrix(0, nrow = ctnum * week, ncol = ctnum)
    for(ct in 1:ctnum){
      A[((ct - 1) * week + 1) : (ct * week), ct] = 1
    }

    # DLM Structure: Choose critical windows for PM2.5 and Temperature
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))
    start.time1 <- sample(1:(Lags - 7), 1)
    start.time2 <- sample(1:(Lags - 7), 1)
    eff1 <- eff2 <- rep(0, Lags)
    eff1[start.time1:(start.time1 + 7)] <- 1
    eff2[start.time2:(start.time2 + 7)] <- 1

    # Data generating process
    # Start a y vector as 0.
    y <- rep(0, n.samp)
    
    # Need to redefine c (renamed to eta1 & eta2) as we have beta1 for binary component and beta2 for negbin component
    beta1 <- rnorm(ncol(data_zi)) # Sample true beta1 from a standard normal
    eta1 <- c(data_zi[, 1:ncol(data_zi)] %*% beta1)     # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          
    eta1_sp <- eta1 + A%*%sp_phi

    phi <- 1 / (1 + exp(-eta1_sp))        # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, phi)           # "At-risk" indicator variable
    nStar <- sum(w)                       # Proportion of at-risk zeros

    # Sample beta2 for negbin component and use the variable w to calculate eta2 (At-risk observations only)
    beta2 <- rnorm(ncol(data_nb))  # Sample true beta2 from a standard normal
    eta2 <- c(data_nb[w == 1, 1:ncol(data_nb)] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Compute f
    truth1Sums <- exposures[[1]][w == 1, ] %*% eff1
    truth2Sums <- exposures[[2]][w == 1, ] %*% eff2
    #f <- truth1Sums + (0.025 * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    int.effect = 0.025
    f <- truth1Sums + (int.effect * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    effect.size <- 1      #/sd(f) # default = 0.1
    f <- f * effect.size

    scale.size = 0.1
    eta2_dlm <- scale.size * (eta2 + f)
    pi <- 1 / (1 + exp(-eta2_dlm))            # Probability of success in negative binomial

    # draw y with eta1, eta2, and f
    r <- 1                                  # Dispersion parameter
    mu <- r * pi / (1 - pi)
    y[w == 1] <- rnbinom(nStar, r, mu = mu)

    # Compute marginal effects
    # Scalar such that the variance of f is 1
    # f <- f * effect.size # Scaling 

    truthInt <- outer(eff1, eff2) * int.effect * effect.size * scale.size
    margDLM1 <- eff1 * effect.size * scale.size + rowSums(truthInt) * mean(exposures[[2]])
    margDLM2 <- colSums(truthInt) * mean(exposures[[1]])

    zeroStr <- length(y[y == 0 & w == 0])/n # y = 0 & classified as structural
    zeroAr <- length(y[y == 0 & w == 1])/n  # y = 0 & classified as at-risk
    nonzero <- length(y[y != 0])/n # y is not zero hence structural
    zeroProp <- length(y[y == 0])/n

    return(list("data" = cbind.data.frame(y, data_zinb),
                "exposures" = exposures,
                "b1" = beta1,
                "b2" = beta2,
                "eff1" = eff1,
                "eff2" = eff2,
                "outer" = outer(eff1, eff2),
                "start.time1" = start.time1, "start.time2" = start.time2,
                "margDLM1" = margDLM1, "margDLM2" = margDLM2,
                "c" = c, "f" = f, "w" = w, "mu" = mu, "phi" = phi,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp,
                "Qinv" = Qinv,
                "eta1" = eta1,
                "areaA" = A, 
                "D" = D, "W" = W,
                "spNodes" = spNodes,
                "spPhi" = sp_phi,
                "AspPhi" = A%*%sp_phi,
                "eta1sp" = eta1_sp))
  }
}
