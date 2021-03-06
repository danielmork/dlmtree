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
                      expList = NULL)
{
  if (!(sim %in% 1:4))
    stop("`sim` must be an integer from 1-4")
  if (is.null(expList)) {
    data(exposureCov)
    cholexp <- t(chol(exposureCov))
    exposureDat <- t(sapply(1:n, function(i) cholexp %*% rnorm(ncol(cholexp))))
    exposureDat <- exposureDat - min(exposureDat)
    expList <- list("e1" = exposureDat[,1:37],
                    "e2" = exposureDat[,38:74],
                    "e3" = exposureDat[,75:111],
                    "e4" = exposureDat[,112:148],
                    "e5" = exposureDat[,149:185])
  }
  Lags <- ncol(expList[[1]])
  n.samp <- min(nrow(expList[[1]]), n)
  idx <- sample(nrow(expList[[1]]), size = n.samp)
  data <- cbind(matrix(rnorm(5 * n.samp), n.samp, 5),
                matrix(rbinom(5 * n.samp, 1, .5), n.samp, 5))
  colnames(data) <- c(paste0("c", 1:5), paste0("b", 1:5))
  params <- rnorm(10)
  c <- c(data[,1:10] %*% params)

  # Sim 1: binary response with single exposure DLM
  if (sim == 1) {
    # center/scale exposure data
    exposures <- lapply(expList, function(i) (i[idx,] - mean(i[idx,])) / sd(i[idx,]))
    # generate random starting time
    start.time1 <- sample(1:(Lags - 7), 1)
    eff1 <- rep(0, Lags)
    eff1[start.time1:(start.time1 + 7)] <- 1
    truth1Sums <- exposures[[1]] %*% eff1
    f <- truth1Sums
    # generate output and set alpha for mean.p
    pp <- (-0.1) * (c + f)
    p <- 1 / (1 + exp(pp - mean(pp) + log(mean.p^(-1) - 1)))
    f <- 0.1 * f
    # draw y
    y <- rbinom(n.samp, 1, p)
    margDLM1 <- eff1 * 0.1
    return(list("dat" = cbind.data.frame(y, data), "params" = params * 0.1,
                "exposures" = exposures, "start.time1" = start.time1,
                "margDLM" = margDLM1, "c" = c, "f" = f, "p" = p))
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
    truth1Sums <- exposures[[1]] %*% eff1
    truth2Sums <- exposures[[2]] %*% eff2
    f <- truth1Sums + (0.025 * truth1Sums * truth2Sums)
    effect.size <- 1 / sd(f)
    f <- f * effect.size
    # draw y
    y <- c + f + rnorm(n.samp, sd = sqrt(error))
    truthInt <- outer(eff1, eff2) * 0.025 * effect.size
    margDLM1 <- eff1 * effect.size + rowSums(truthInt) * mean(exposures[[2]])
    margDLM2 <- colSums(truthInt) * mean(exposures[[1]])
    return(list("dat" = cbind.data.frame(y, data),
                "exposures" = exposures, "params" = params,
                "start.time1" = start.time1, "start.time2" = start.time2,
                "margDLM1" = margDLM1, "margDLM2" = margDLM2,
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


  # Sim 4: continuosu response to test exposure selection using one exposure
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
