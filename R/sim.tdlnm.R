#' sim.tdlnm
#' 
#' @description Simulation scenarios to accompany TDLNM
#'
#' @param effect character (A - D) indicating simulation scenario
#' @param error.to.signal scalar value setting error: sigma^2/var(f)
#'
#' @return
#' @export
#'
sim.tdlnm <- function(effect = "A", error.to.signal = 1)
{
  data("pm25Exposures")
  pm25Exposures <- log(pm25Exposures[which(pm25Exposures$S == "Colorado"),-c(1:2)])[,1:37]
  n.samp <- nrow(pm25Exposures)

  data <- cbind(matrix(rnorm(5*n.samp), n.samp, 5),
                matrix(rbinom(5*n.samp, 1, .5), n.samp, 5))
  colnames(data) <- c(paste0("c", 1:5), paste0("b", 1:5))
  params <- rnorm(10)
  c <- c(data[,1:10] %*% params)

  # Piecewise constant effect ----
  if (effect == "A") {
    cenval = 1
    dlnm.fun <- function(exposure.data, cenval, sum = T) {
      if (sum) {
        -rowSums(exposure.data[, 11:15] > 2)
      } else {
        dlnm <- t(sapply(1:nrow(exposure.data), function(i) {
          sapply(1:ncol(exposure.data), function(j) {
            ifelse(j %in% 11:15, ifelse(exposure.data[i, j] > 2, -1, 0), 0)
          })}))
        colnames(dlnm) <- paste0("Time", 1:ncol(exposure.data))
        return(dlnm)
      }
    }

  # Linear effect ----
  } else if (effect == "B") {
    cenval = 1
    dlnm.fun <- function(exposure.data, cenval, sum = T) {
      if (sum) {
        rowSums((cenval - exposure.data[, 11:15]))
      } else {
        dlnm <- t(sapply(1:nrow(exposure.data), function(i) {
          sapply(1:ncol(exposure.data), function(j) {
            ifelse(j %in% 11:15, (cenval - exposure.data[i, j]), 0)
          })}))
        colnames(dlnm) <- paste0("Time", 1:ncol(exposure.data))
        return(dlnm)
      }
    }

  # Logistic effect, piecewise in time ----
  } else if (effect == "C") {
    cenval = 1
    dlnm.fun <- function(exposure.data, cenval, sum = T) {
      flogistic <- function(x) ((1/(1+exp(5*(x-2.5)))) - 1)

      if (sum) {
        rowSums(flogistic(exposure.data[, 11:15]))
      } else {
        dlnm <- t(sapply(1:nrow(exposure.data), function(i) {
          sapply(1:ncol(exposure.data), function(j) {
            ifelse(j %in% 11:15, flogistic(exposure.data[i, j]), 0)
          })}))
        colnames(dlnm) <- paste0("Time", 1:ncol(exposure.data))
        return(dlnm)
      }
    }

  # Logistic effect - smooth in time ----
  } else if (effect == "D") {
    cenval = 1
    dlnm.fun <- function(exposure.data, cenval, sum = T) {
      flogistic <- function(x) ((1/(1+exp(5*(x-2.5)))) - 1)

      ftime <- function(t) (exp(-.0025 * (t-13)^4))

      f <- function(x, t) (flogistic(x) * ftime(t))

      if (sum) {
        sapply(1:nrow(exposure.data), function(i) {
          sum(f(exposure.data[i,,drop=T], (1:ncol(exposure.data))))
        })
      } else {
        dlnm <- t(sapply(1:nrow(exposure.data), function(i) {
          sapply(1:ncol(exposure.data), function(j) {
            f(exposure.data[i, j], j)
          })}))
        colnames(dlnm) <- paste0("Time", 1:ncol(exposure.data))
        return(dlnm)
      }
    }
  }

  f <- dlnm.fun(pm25Exposures, cenval, T)
  e <- rnorm(length(f), 0, sqrt(var(f) * error.to.signal))
  y <- c + f + e


  return(list("dat" = cbind.data.frame(y, data),
              "exposure" = as.matrix(pm25Exposures),
              "dlnm.fun" = dlnm.fun,
              "cenval" = cenval,
              "params" = params,
              "c" = c,
              "f" = f))
}
