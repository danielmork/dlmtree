#' sim.hdlmm
#'
#' @title Creates simulated data for HDLM & HDLMM
#' @description Method for creating simulated data for HDLM & HDLMM
#'
#' @param sim character (A - E) specifying simulation scenario
#' @param n sample size
#' @param error positive scalar specifying error variance for Gaussian response
#' @param effect.size the effect size of the window of susceptibility
#' @param exposure.data exposure data. A matrix of exposure data for simulation A, B, C and a named list of exposure data for simulation D, E
#'
#' @details Simulation scenarios:
#' - Scenario A: Two subgroups with early/late windows determined by continuous and binary modifiers 
#' - Scenario B: Two subgroups with scaled effect determined by a continuous modifier
#' - Scenario C: No heterogeneity i.e., same effect on all individuals
#' - Scenario D: Three subgroups with three corresponding exposures. Subgroups are determined by continuous and binary modifiers
#' - Scenario E: Two subgroups with two exposures. First group is associated with the scaled main effect and lagged interaction 
#' while the second group is only associated with the scaled main effect, no interaction.
#' @md
#'
#' @examples
#' sim.hdlmm(sim = "A", n = 1000)
#'
#' @returns Simulated data and true parameters
#' @export
#'
sim.hdlmm <- function(sim = "A",
                      n = 1000,
                      error = 1,
                      effect.size = 1,
                      exposure.data = NULL)
{
  # Exposure data availability check
  if (is.null(exposure.data)) {
    if (sim %in% c("A", "B", "C")) {
      data("pm25Exposures", envir = environment())
      exposure.data <- sapply(3:39, function(i) pm25Exposures[,i])
    } else {
      data("coExp", envir = environment())
      if (sim == "D") {
        exposure.data <- list("e1" = coExp[,1:37],
                              "e2" = coExp[,38:74],
                              "e3" = coExp[,75:111])
      } else {
        exposure.data <- list("e1" = coExp[,1:37],
                              "e2" = coExp[,38:74])
      }
    } 
  }

  # Exposure effect 
  if (sim %in% c("A", "B", "C")) {
    exposure.data <- (exposure.data - mean(exposure.data)) / sd(exposure.data)
    pX            <- ncol(exposure.data)
    n.samp        <- min(nrow(exposure.data), n)
    exposure.data <- exposure.data[sample(nrow(exposure.data), n.samp),]
  } else { # Multiple exposures
    pX            <- ncol(exposure.data[[1]])
    n.samp        <- min(nrow(exposure.data[[1]]), n)
    idx           <- sample(nrow(exposure.data[[1]]), size = n.samp)
    exposure.data <- lapply(exposure.data, function(i) (i[idx, ] - mean(i[idx, ])) / sd(i[idx, ]))
  } 

  # Fixed effect
  dat <- data.frame(rnorm(n.samp), rbinom(n.samp, 1, 0.5), runif (n.samp),
                    matrix(rnorm(n.samp * 5), n.samp, 5),
                    matrix(rbinom(n.samp * 5, 1, 0.5), n.samp, 5))
  colnames(dat) <- c("mod_num", "mod_bin", "mod_scale", paste0("c", 1:5), paste0("b", 1:5))

  sd.f   <- 1
  dlmFun <- function(dat.row) {}
  if (sim == "A") {
    dlmFun <- function(dat.row) {
      if (dat.row$mod_num > 0) {  
        if (dat.row$mod_bin == 1) {
          c(rep(0, 10), rep(1, 8), rep(0, pX - 18)) / sd.f 
        } else {
          c(rep(0, 16), rep(1, 8), rep(0, pX - 24)) / sd.f 
        }
      } else {
        rep(0, pX)
      }
    }

    fixedIdx <- list(which(dat$mod_num > 0 & dat$mod_bin == 1),   
                     which(dat$mod_num > 0 & dat$mod_bin == 0),  
                     which(dat$mod_num < 0))  
  } else if (sim == "B") {
    dlmFun <- function(dat.row) {
      if (dat.row$mod_num > 0) {
        c(rep(0, 10), rep(dat.row$mod_scale, 8), rep(0, pX - 18)) / sd.f
      } else {
        rep(0, pX)
      }
    }
    fixedIdx <- list(which(dat$mod_num > 0 & dat$mod_scale < 0.25),
                     which(dat$mod_num > 0 & dat$mod_scale > 0.25 & dat$mod_scale < 0.5),
                     which(dat$mod_num > 0 & dat$mod_scale > 0.5 & dat$mod_scale < 0.75),
                     which(dat$mod_num > 0 & dat$mod_scale > 0.75),
                     which(dat$mod_num < 0))
  } else if (sim == "C") {
    start.time <- sample(1:(pX - 7), 1)
    dlmFun <- function(dat.row) {
      d <- rep(0, pX)
      d[start.time:(start.time + 7)] <- 1
      d
    }
    fixedIdx <- list(1:nrow(dat))
  } else if (sim == "D") {  
    start.time1 <- sample(1:(pX - 7), 1)
    start.time2 <- sample(1:(pX - 7), 1)
    start.time3 <- sample(1:(pX - 7), 1)
    dlmFun <- function(dat.row) {
      if (dat.row$mod_num > 0) { 
        if (dat.row$mod_bin == 1) {
          e1 <- rep(0, pX)
          e1[start.time1:(start.time1 + 7)] <- effect.size
          e1
        } else {
          e2 <- rep(0, pX)
          e2[start.time2:(start.time2 + 7)] <- effect.size
          e2
        }
      } else {
        e3 <- rep(0, pX)
        e3[start.time3:(start.time3 + 7)] <- effect.size
        e3
      }
    }
    fixedIdx <- list(which(dat$mod_num > 0 & dat$mod_bin == 1),   
                     which(dat$mod_num > 0 & dat$mod_bin == 0),  
                     which(dat$mod_num < 0))               

  } else if (sim == "E") {
    start.time11  <- sample(1:(pX - 7), 1)  # Subgroup 1 - main t4
    start.time12  <- sample(1:(pX - 7), 1)  # Subgroup 1 - interaction t6
    start.time2   <- sample(1:(pX - 7), 1)  # Subgroup 2 - main t5
    
    eff11 <- eff12 <- eff2 <- rep(0, pX)
    eff11[start.time11:(start.time11 + 7)]  <- effect.size
    eff12[start.time12:(start.time12 + 7)]  <- effect.size
    eff2[start.time2:(start.time2 + 7)]     <- effect.size
    int.size <- 0.025

    dlmFun <- function(dat.row) {
      if (dat.row$mod_num > 0) {  
        list("e1" = eff11, "e2" = eff12) 
      } else {
        eff2 
      }
    }
    
    fixedIdx <- list(which(dat$mod_num > 0), 
                     which(dat$mod_num <= 0)) 
  } 

  # Calculate the DLM effect, f
  f <- rep(NA, n.samp)

  if (!(sim %in% c("D", "E"))) {
    f <- sapply(1:n.samp, function(i) sum(exposure.data[i,] * dlmFun(dat[i, , drop = FALSE])))

  } else if (sim == "D") {
    for (group in 1:length(fixedIdx)) { 
      indices <- fixedIdx[[group]]
      n.group <- length(indices)
      
      for (i in 1:n.group) {
        currentIdx    <- indices[i]
        f[currentIdx] <- sum(exposure.data[[group]][currentIdx, ] * dlmFun(dat[currentIdx, , drop = FALSE]))
      } 
    }
  } else if (sim == "E") {
    for (group in 1:length(fixedIdx)) { 
      indices <- fixedIdx[[group]]
      n.group <- length(indices)
      
      for (i in 1:n.group) {
        currentIdx <- indices[i]

        if (group == 1) {
          effList   <- dlmFun(dat[currentIdx, , drop = FALSE])

          e1_effect <- sum(exposure.data[[1]][currentIdx, ] * effList$e1)
          e2_effect <- sum(exposure.data[[2]][currentIdx, ] * effList$e2)

          f[currentIdx] <- dat[currentIdx, , drop = FALSE]$mod_scale * e1_effect + int.size * (e1_effect * e2_effect) # z1e1 + e1xe2 
        } else {
          f[currentIdx] <- dat[currentIdx, , drop = FALSE]$mod_scale * sum(exposure.data[[1]][currentIdx, ] * dlmFun(dat[currentIdx, , drop = FALSE])) # z1e1
        } 
      }
    }
  } 

  # Scale f
  sd.f  <- sd(f)
  f     <- f / sd.f

  # y sample
  params  <- rnorm(13)
  c       <- as.matrix(dat) %*% params
  dat$y   <- c + f + rnorm(n.samp, sd = sqrt(error))

  return(list("dat" = dat, 
              "exposures" = exposure.data,
              "f" = f, 
              "c" = c,
              "params" = params, 
              "sd.f" = sd.f,
              "fixedIdx" = fixedIdx,
              "dlmFun" = dlmFun))
}
