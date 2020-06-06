#' dlnmMod
#'
#' @param X Exposure data
#' @param SE Standard error of exposure data, can be left null
#' @param X.bins Interger value specifying the number of evenly spaced quantiles
#'  to split exposures into (0 indicates to use all values as possible splits)
#'  or list("type" = c("values", quantiles"), split.vals = numeric vector)
#'
#' @description Modifiers environment for DLNM
#'
#' @return Environment with objects and functions
#' @importFrom stats quantile rbinom
#'
#' @examples
dlnmMod <- function(X, SE = NULL, X.bins)
{
  X <- as.matrix(X)
  if (!is.null(SE)) {
    SE <- as.matrix(SE)
  }

  # X.bins is specified as a list of arguments, can be split at
  # quantiles or specific values, or equally spaced quantiles or values
  if (!is.list(X.bins)) {
    if (length(X.bins) > 1)
      Xsplits <- X.bins
    else if (is.numeric(X.bins))
      Xsplits <- unique(sort(quantile(c(X), 1:(X.bins-1)/X.bins)))
    else if (X.bins == Inf)
      Xsplits <- unique(sort(c(X)))
    else
      Xsplits <- c()

  } else {
    if (X.bins$type == "values") {
      Xsplits <- unique(sort(X.bins$split.vals))
      Xsplits <- Xsplits[which(Xsplits > min(X) & Xsplits < max(X))]

    } else if (X.bins$type == "quantiles") {
      if (any(X.bins$split.vals < 0) | any(X.bins$split.vals > 1))
        stop("When using quantiles, split.vals must be between 0 and 1")
      Xsplits <- unique(sort(quantile(c(X), X.bins$split.vals)))
    }
  }

  if (length(Xsplits) == 0) {
    Xscale <- diff(range(X))
    Xmean <- mean(X)
    X <- (X - Xmean) / Xscale
  } else {
    Xscale <- sqrt(nrow(X)) * ncol(X)
    Xmean <- 0
  }


  # Pre-count exposure splits for all times (these happen often)
  Xcounts <- list()
  Tsums <- list()
  ZtX <- list()
  VgZtX <- list()
  preset.counts <- function(Z, Vg) {
    p.counts <- rep(1 / sqrt(nrow(X)), nrow(X))
    remove <- c()

    if (length(Xsplits) == 0) {
      Xcounts[[1]] <<- rowSums(X)
      for (i in 1:(ncol(X)-1)) {
        if (i == 1) {
          Tsums[[i]] <<- list(X[,1], rowSums(X[,2:ncol(X)]))
        } else if (i == ncol(X)-1) {
          Tsums[[i]] <<- list(rowSums(X[,1:(ncol(X)-1)]), X[,ncol(X)])
        } else {
          Tsums[[i]] <<- list(rowSums(X[,1:i]), rowSums(X[,(i+1):ncol(X)]))
        }
        ZtX[[i]] <<- list(crossprod(Z, Tsums[[i]][[1]]), crossprod(Z, Tsums[[i]][[2]]))
        VgZtX[[i]] <<- list(crossprod(Vg, ZtX[[i]][[1]]), crossprod(Vg, ZtX[[i]][[2]]))
      }
      return(invisible())
    }

    for (i in 1:length(Xsplits)) {

      if (is.null(SE))
        nc <- nodeCount(X, Z, Vg, p.counts, -Inf, Xsplits[i], 1, ncol(X))
      else
        nc <- nodeCountSE(X, SE, Z, Vg, p.counts, -Inf, Xsplits[i], 1, ncol(X))

      if (!nc$Empty) {
        Xcounts[[length(Xcounts) + 1]] <<- nc$Count
        ZtX[[length(Xcounts)]] <<- nc$ZtX
        VgZtX[[length(Xcounts)]] <<- nc$VgZtX
      } else {
        remove <- c(remove, i)
      }
    }

    if (length(remove) > 0)
      Xsplits <<- Xsplits[-remove]

    Xcounts[[length(Xcounts) + 1]] <<- rep(1 / sqrt(nrow(X)), nrow(X))
    ZtX[[length(Xcounts)]] <<- crossprod(Z, Xcounts[[length(Xcounts)]])
    VgZtX[[length(Xcounts)]] <<- crossprod(Vg, ZtX[[length(Xcounts)]])
  }


  Tsplits <- 1:(ncol(X) - 1)
  p.x <- ifelse(length(Xsplits) > 0, 1 / length(Xsplits), 0)
  p.t <- 1 / length(Tsplits)

  # Count available splits remaining
  avail.mod = function(node) {
    return(list("Xs" = Xsplits[which(Xsplits > node$xmin & Xsplits < node$xmax)],
                "Ts" = Tsplits[which(Tsplits >= node$tmin & Tsplits < node$tmax)]))
  }

  # Calculate log probability of rule
  log.pRule = function(node) {
    if (is.null(node$avail.mod))
      node$avail.mod <- avail.mod(node)

    lX <- length(node$avail.mod$Xs)
    lT <- length(node$avail.mod$Ts)

    if (!is.null(node$xsplit))
      return(log(p.x / (lX * p.x + lT * p.t)))
    else
      return(log(p.t / (lX * p.x + lT * p.t)))

  }

  # Propose new split for current node
  propose.split = function(node) {
    lX <- length(node$avail.mod$Xs)
    lT <- length(node$avail.mod$Ts)

    if (is.null(node$avail.mod))
      node$avail.mod <- avail.mod(node)

    if (lX == 0 && lT == 0) {
      return(list("Success" = FALSE))

    } else if (lX == 0 | lT == 0) {
      if (lT > 0)
        m <- "T"
      else
        m <- "X"

    } else {
      if (rbinom(1, 1, lX * p.x / (lX * p.x + lT * p.t)))
        m <- "X"
      else
        m <- "T"
    }

    if (m == "X") {
      return(list("Success" = TRUE,
                  "logprule" = log(p.x / (lX * p.x + lT * p.t)),
                  "mod" = m,
                  "split" = node$avail.mod$Xs[sample.int(length(node$avail.mod$Xs), 1)]))
    } else {
      return(list("Success" = TRUE,
                  "logprule" = log(p.t / (lX * p.x + lT * p.t)),
                  "mod" = m,
                  "split" = node$avail.mod$Ts[sample.int(length(node$avail.mod$Ts), 1)]))
    }
  }

  environment()
}
