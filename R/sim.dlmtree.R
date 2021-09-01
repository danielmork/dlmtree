sim.dlmtree <- function(sim = "A",
                        error = 1,
                        n = 1000,
                        exposure.data = NULL)
{
  if (is.null(exposure.data)) {
    data("pm25Exposures")
    exposure.data = sapply(3:39, function(i) pm25Exposures[,i])
  }

  exposure.data <- (exposure.data - mean(exposure.data)) / sd(exposure.data)
  pX <- ncol(exposure.data)
  n.samp <- min(nrow(exposure.data), n)
  exposure.data <- exposure.data[sample(nrow(exposure.data), n.samp),]
  dat <- data.frame(rnorm(n.samp), rbinom(n.samp, 1, 0.5), runif(n.samp),
                    matrix(rnorm(n.samp * 5), n.samp, 5),
                    matrix(rbinom(n.samp * 5, 1, 0.5), n.samp, 5))
  colnames(dat) <- c("mod_num", "mod_bin", "mod_scale",
                     paste0("c", 1:5), paste0("b", 1:5))

  sd.f <- 1
  dlmFun <- function() {}
  if (sim == "A") {
    dlmFun <- function(dat.row) {
      if (dat.row$mod_num > 0) {
        if (dat.row$mod_bin == 1)
          c(rep(0, 10), rep(1, 8), rep(0, pX - 18)) / sd.f
        else
          c(rep(0, 16), rep(1, 8), rep(0, pX - 24)) / sd.f
      } else
        rep(0, pX)
    }
    fixedIdx <- list(which(dat$mod_num > 0 & dat$mod_bin == 1),
                     which(dat$mod_num > 0 & dat$mod_bin == 0),
                     which(dat$mod_num < 0))
  } else if (sim == "B") {
    dlmFun <- function(dat.row) {
      if (dat.row$mod_num > 0) {
        c(rep(0, 10), rep(dat.row$mod_scale, 8), rep(0, pX - 18)) / sd.f
      } else
        rep(0, pX)
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
  }

  f <- sapply(1:n.samp, function(i) sum(exposure.data[i,] * dlmFun(dat[i,,drop=F])))
  sd.f <- sd(f)
  f <- f / sd.f

  params <- rnorm(13)
  c <- as.matrix(dat) %*% params

  dat$y <- c + f + rnorm(n.samp, sd = sqrt(error))

  return(list("dat" = dat, "exposure.dat" = exposure.data,
              "f" = f, "c" = c,
              "params" = params, "sd.f" = sd.f,
              "fixedIdx" = fixedIdx,
              "dlmFun" = dlmFun))
}
