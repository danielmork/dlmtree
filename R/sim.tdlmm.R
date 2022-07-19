#' sim.tdlmm
#'
#' @description Simulation scenarios to accompany TDLM/TDLMM
#'
#' @param sim integer (1-4) specifying simulation scenario
#' @param error positive scalar specifying error variance
#' @param mean.p scalar between zero and one specifying mean probability
#' for simulation scenario one
#' @param n.exp number of exposures for simulation scenarios three and four
#' @param prop.active proportion of active exposures for simulation scenario
#' three
#' @param n sample size for simulation
#' @param expList named list of exposure data
#'
#' @return
#' @export
#'
sim.tdlmm <- function(sim = 1,
                      error = 10,
                      mean.p = 0.5,
                      n.exp = 25,
                      prop.active = 0.05,
                      n = 5000,
                      spatial = FALSE,
                      expList = NULL,
                      areaN = NULL,
                      rho = NULL,
                      spTau = NULL)
{
  if (!(sim %in% 1:7))
    stop("`sim` must be an integer from 1-7")

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
                matrix(rbinom(5 * n.samp, 1, .5), n.samp, 5))   # 5 binary predictors (n.samp x 5)                                                               # => data: (n.samp x 10)
  colnames(data) <- c(paste0("c", 1:5), paste0("b", 1:5))       # Name the columns with c1 - c5 and b1 - b5
                                                                # c for continuous, b for binary
  params <- rnorm(10)                                           # Sample true beta from a standard normal
  c <- c(data[,1:10] %*% params)                                # Compute c: zT * gamma (nx10) x (10x1) = (nx1)

  # Sample spatial geoid
  # Currently using a fake spatial geoid from Denver county
  # if(spatial){
  #   fakeGeoid = c("08031003603", "08031004103", "08031004104", "08031003701", "08031004201", "08031004202",
  #                  "08031003702", "08031003703", "08031003300", "08031004301", "08031004302", "08031004304")
  #   geoid = sample(fakeGeoid, n, replace = T)
  # } else {
  #   geoid = NULL
  # }


  if(spatial){
    # If the parameter, spatial == true, create a fake adjacency matrix
    Q = diag(rep(1, areaN)) # Identity matrix
    A = matrix(0, nrow = n, ncol = areaN)


    if(n < areaN){
      stop("Sample size cannot be less than the total number of areal units")
    }

    #print("entering spatial")
    allAdj = FALSE # Flag to check if all of the observations has at least one neighbor

    #print("creating lower triangle")
    while(!allAdj){
      # Randomly choose edges # between 1 ~ n/2 and create W
      nConnect = sample(round(areaN/2):round(areaN*(areaN-1)/2), 1)
      W = lowerAdj(areaN, nConnect) 
      W = mirrorMat(W)

      # If all of the observations has at least one neighbor
      if(all(rowSums(W) > 0)){ 
        allAdj = TRUE
      }
    }
    
    #print("Computing D, Q")
    # Now that we have created a lower matrix of adjacency, generate W, D, Q
    D <- diag(rowSums(W))

    # ICAR var-cov
    Qsim <- spTau * (D - rho*W)
    Qsiminv <- chol2inv(chol(Qsim))

    # Compute the spatial edge nodes
    spNodes = spatial_pairwise(W)
    
    # Sorting matrix of 0 & 1: A x phi = (n x areaN) (areaN x 1) -> areaIndx of each observation is the column of A
    # Assign the areas to the observations
    for(i in 1:n){
      mod = i %% areaN
      A[i,  mod + 1] = 1 # There is no index for mod 0 so add one
    }
  }

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

  # Sim 5: Zero-inflated negative binomial with single exposure DLM (temporarily here as "p" is not defined)
  if (sim == 5) {
    data <- cbind(matrix(rnorm(5 * n.samp), n.samp, 5),         # 5 guassian predictors (n.samp x 2)
                matrix(rbinom(5 * n.samp, 1, .5), n.samp, 5))   # 5 binary predictors (n.samp x 2)                                                               # => data: (n.samp x 10)
    colnames(data) <- c(paste0("c", 1:5), paste0("b", 1:5))     # Name the columns with c1 - c2 and b1 - b2
                                                                # c for continuous, b for binary
    #params <- rnorm(10)                                          # Sample true beta from a standard normal
    #c <- c(data[,1:10] %*% params)                               # Compute c: zT * gamma (nx4) x (4x1) = (nx1)

    # center/scale exposure data
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))
    # generate random starting time
    start.time1 <- sample(1:(Lags - 7), 1)      # sample from 1 ~ 30
    eff1 <- rep(0, Lags)                        # eff1 = a vector of 37 zeros
    eff1[start.time1:(start.time1 + 7)] <- 1    # Set 7 weeks of eff1 as 1

    # Start a y vector as 0.
    y <- rep(0, n.samp)

    # Need to redefine c (renamed to eta1 & eta2) as we have beta1 for binary component and beta2 for negbin component
    beta1 <- rnorm(10) # Sample true beta1 from a standard normal
    eta1 <- c(data[, 1:10] %*% beta1)     # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          

    phi <- 1 / (1 + exp(-eta1))           # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, phi)           # "At-risk" indicator variable
    nStar <- sum(w)                       # Proportion of structural zeros

    # Sample beta2 for negbin component and use the variable w to calculate eta2 (At-risk observations only)
    beta2 <- rnorm(10)  # Sample true beta2 from a standard normal
    eta2 <- c(data[w == 1, 1:10] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Update f to use only at-risk
    truth1Sums <- exposures[[1]][w == 1, ] %*% eff1       # A single exposure (n* x 37) %*% eff1 (37 x 1) = (n* x 1)
    f <- truth1Sums # (n x 1)
    effect.size = 1 # 1/sd(f)
    f <- f * effect.size

    scale.size = 0.25
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
    return(list("dat" = cbind.data.frame(y, data), 
                "params" = params * 0.25, # return the true coefficient of 10 parameters
                "exposures" = exposures, 
                "start.time1" = start.time1,
                "margDLM1" = margDLM1, 
                "eta1" = eta1, "eta2" = eta2,
                "beta1" = beta1, "beta2" = beta2,
                "phi" = phi, # At-risk probability
                "mu" = mu, "pi" = pi, "r" = r, "f" = f, "w" = w,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp))
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

  # Sim 6: ZINB with 5 exposures with/without interaction
  if(sim == 6){
    # DLM Structure: Choose critical windows for PM2.5 and NO2
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
    beta1 <- rnorm(10) # Sample true beta1 from a standard normal
    eta1 <- c(data[, 1:10] %*% beta1)     # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)          

    phi <- 1 / (1 + exp(-eta1))           # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, phi)           # "At-risk" indicator variable
    nStar <- sum(w)                       # Proportion of at-risk zeros

    # Sample beta2 for negbin component and use the variable w to calculate eta2 (At-risk observations only)
    beta2 <- rnorm(10)  # Sample true beta2 from a standard normal
    eta2 <- c(data[w == 1, 1:10] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Compute f
    truth1Sums <- exposures[[1]][w == 1, ] %*% eff1
    truth2Sums <- exposures[[2]][w == 1, ] %*% eff2
    #f <- truth1Sums + (0.025 * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    int.effect = 0.025
    f <- truth1Sums + (int.effect * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    effect.size <- 1
    f <- f * effect.size

    scale.size = 0.25
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

    return(list("dat" = cbind.data.frame(y, data),
                "exposures" = exposures, "params" = params,
                "eff1" = eff1,
                "eff2" = eff2,
                "outer" = outer(eff1, eff2),
                "truthInt" = truthInt,
                "start.time1" = start.time1, "start.time2" = start.time2,
                "margDLM1" = margDLM1, "margDLM2" = margDLM2,
                "c" = c, "f" = f, "w" = w, "mu" = mu, "phi" = phi, "eta1" = eta1,
                "b1" = beta1,
                "b2" = beta2,
                "zeroStr" = zeroStr,
                "zeroAr" = zeroAr,
                "nonzero" = nonzero,
                "zeroProportion" = zeroProp))
  }

  # Sim 7: ZINB + TDLMM + Spatial random effect
  if(sim == 7){
    # # Spatial effect simulation
    # stateCode <- substr(geoid[1], 1, 2)          # Could be a vector but usually one number
    # countiesCode <- unique(substr(geoid, 3, 5))  # A vector

    # # With geoid, create a SpatialPolygonsDataFrame (SPDF) object
    # ctspatial <- tigris::tracts(state = stateCode, county = countiesCode, refresh = T, class = "sp") 

    # # geoid_df & sp_data
    # geoid_df = data.frame("geoid" = geoid)
    # geoid_df$state = substr(geoid_df$geoid, 1, 2)  
    # geoid_df$county = substr(geoid_df$geoid, 3, 5)  
    # geoid_df$tract = substr(geoid_df$geoid, 6, 11)  

    # sp_data = cbind(data, geoid_df)

    # # Merge the spatial information with the cityfile
    # spatial_data <- sp::merge(ctspatial, sp_data, by.y = "geoid", by.x = "GEOID", all.x = F, all.y = T, duplicateGeoms = TRUE)

    # # Get an adjacency matrix
    # adj <- poly2nb(spatial_data)

    # #adj <- poly2nb(spatial_data, row.names = spatial_data$GEOID)
    # coord_df = coordinates(spatial_data)
    # row.names(coord_df) = spatial_data$GEOID

    # W <- nb2mat(adj, style = "B") # Adjacency matrix
    # D <- diag(rowSums(W))         # Diagonal matrix

    # Compute variance-covariance matrix for phi (CAR model)
    #alpha = 0.5 # alpha = 0 implies no spatial effect. Note that alpha = 1 makes Q non-invertible
    #Q <- chol2inv(chol(D - alpha*W))
    spatial_phi = c(mvrnorm(n = 1, rep(0, areaN), Qsiminv)) # Simulate spatial effect -> need MASS package

    # DLM Structure: Choose critical windows for PM2.5 and NO2
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
    beta1 <- rnorm(10) # Sample true beta1 from a standard normal
    eta1 <- c(data[, 1:10] %*% beta1)     # Compute eta1: xT * beta1 (n x 10)x(10 x 1) = (nx1)   
    eta1_sp <- eta1 + A%*%spatial_phi
    #eta1_sp <- A%*%spatial_phi
    spRatio <- eta1 / A%*%spatial_phi

    phi <- 1 / (1 + exp(-eta1_sp))        # 1 - P(structural zero)
    w <- rbinom(n.samp, 1, phi)           # "At-risk" indicator variable
    nStar <- sum(w)                       # Proportion of at-risk zeros

    # Sample beta2 for negbin component and use the variable w to calculate eta2 (At-risk observations only)
    beta2 <- rnorm(10)  # Sample true beta2 from a standard normal
    eta2 <- c(data[w == 1, 1:10] %*% beta2)     # Compute eta2: xT * beta2 (n x 10)x(10 x 1) = (nStarx1)

    # Compute f
    truth1Sums <- exposures[[1]][w == 1, ] %*% eff1
    truth2Sums <- exposures[[2]][w == 1, ] %*% eff2
    #f <- truth1Sums + (0.025 * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    int.effect = 0.025
    f <- truth1Sums + (int.effect * truth1Sums * truth2Sums) # PM2.5 main effect + interaction effect
    effect.size <- 1      #/sd(f) # default = 0.1
    f <- f * effect.size

    scale.size = 0.25
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

    return(list("dat" = cbind.data.frame(y, data),
                "exposures" = exposures, "params" = params,
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
                #"geoid" = geoid, 
                "Qinv" = Qsiminv,
                "spRatio" = spRatio,
                "eta1" = eta1,
                "areaA" = A, 
                "D" = D, "W" = W,
                "spNodes" = spNodes,
                "spN" = areaN,
                "spPhi" = spatial_phi,
                "spPhiN" = A%*%spatial_phi,
                "eta1sp" = eta1_sp))
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


}
