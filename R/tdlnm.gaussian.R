#' tdlnm.gaussian
#'
#' @param model Environment with model data and control parameters
#'
#' @description Sets up and runs MCMC for tdlnm model with Gaussian errors.
#' Modifies environment, no return.
#'
#' @return invisible()
#'
#' @examples
tdlnm.gaussian <- function(model)
{
  ctr <- new.env(hash = T)
  ctr$n <- length(model$Y)
  ctr$pZ <- ncol(model$Z)
  ctr$pX <- ncol(model$Mo$X)

  # + Control covariance matrices ----
  ctr$Vg <- solve(crossprod(model$Z) + diag(ctr$pZ) / 100000)
  ctr$VgChol <- chol(ctr$Vg)
  model$Mo$preset.counts(model$Z, ctr$Vg)
  if (length(model$Mo$Xsplits) > 0) {
    ctr$X1 <- rep(1/sqrt(ctr$n), ctr$n)
  } else {
    ctr$X1 <- rowSums(model$Mo$X)
  }
  ctr$XX1 <- as.double(crossprod(ctr$X1))
  ctr$ZtX1 <- crossprod(model$Z, ctr$X1)
  ctr$VgZtX1 <- crossprod(ctr$Vg, ctr$ZtX1)

  # + Initial estimates ----
  B0 <- controlEst(model$Y, model$Z, ctr$Vg, ctr$VgChol, 1, 0, 0)
  ctr$sigma2 <- B0$sigma2
  ctr$xi.sigma2 <- B0$xi.sigma2
  ctr$xi.omega <- 2
  ctr$omega <- 1/rgamma(1, (model$ctr$n.trees + 1) / 2, 1 / ctr$xi.omega)
  ctr$xi.tau <- rep(2, model$ctr$n.trees)
  ctr$tau <- 1/rgamma(model$ctr$n.trees, 1, 1 / ctr$xi.tau)
  ctr$TLiT <- rep(0, model$ctr$n.trees)
  ctr$n.term <- rep(1, model$ctr$n.trees)
  ctr$fhat <- rep(0, ctr$n)
  ctr$Rmat <- lapply(1:(model$ctr$n.trees), function(i) rep(0, ctr$n))
  ctr$R <- ctr$Rmat[[1]]

  # + Create trees ----
  new.tree <- dlnmNode(xmin = -Inf, xmax = Inf, tmin = 1, tmax = ncol(model$Mo$X))
  new.tree$update.subnodes(model$Mo); new.tree$counts <- ctr$X1;
  new.tree$ZtX <- ctr$ZtX1; new.tree$VgZtX <- ctr$VgZtX1
  new.tree$update = F
  ctr$trees <- lapply(1:model$ctr$n.trees, function(i) new.tree$copy())

  # + MCMC Control ----
  lrunif <- log(runif(model$ctr$n.burn + model$ctr$n.iter))

  # + Logs ----
  dgn <- new.env(hash = T)
  dgn$DLM <- lapply(1:floor(model$ctr$n.iter/model$ctr$n.thin), function(i) matrix())
  dgn$gamma <- matrix(0.0, floor(model$ctr$n.iter/model$ctr$n.thin), ctr$pZ)
  dgn$sigma2 <- rep(0.0, floor(model$ctr$n.iter/model$ctr$n.thin))
  dgn$fhat <- rep(0.0, ctr$n)
  if (model$ctr$diagnostics) {
    dgn$omega <- rep(0.0, floor(model$ctr$n.iter/model$ctr$n.thin))
    dgn$tau <- matrix(0.0, floor(model$ctr$n.iter/model$ctr$n.thin), model$ctr$n.trees)
    # dgn$alpha <- rep(0.0, floor(model$ctr$n.iter/model$ctr$n.thin))
    dgn$accept <- lapply(1:(model$ctr$n.trees*floor(model$ctr$n.iter/model$ctr$n.thin)),
                         function(i) list())
    dgn$term.nodes <- matrix(0L, floor(model$ctr$n.iter/model$ctr$n.thin), model$ctr$n.trees)
  }


  ############################################################################ #
  #### MCMC ####
  ############################################################################ #
  if (model$ctr$verbose & !model$ctr$par) {
    cat("Running model '.' = 100 iterations\n")
  }
  ctr$time <- proc.time()[3]
  for (b in 1:(model$ctr$n.burn + model$ctr$n.iter)) {
    ctr$b <- b
    if (ctr$b > model$ctr$n.burn &
        (ctr$b - model$ctr$n.burn) %% model$ctr$n.thin == 0)
      ctr$record <- (ctr$b - model$ctr$n.burn)/model$ctr$n.thin
    else
      ctr$record <- 0

    # + Update all trees ----
    dlnmTreeMCMC(model, ctr, dgn)

    # + Update control estimates ----
    ctr$fhat <- Reduce("+", ctr$Rmat)
    ctr$R <- model$Y - ctr$fhat
    B0 <- controlEst(ctr$R, model$Z, ctr$Vg, ctr$VgChol,
                     ctr$sigma2, ctr$sum.term, ctr$sum.f.exp / ctr$omega)
    ctr$sigma2 <- B0$sigma2
    ctr$xi.sigma2 <- B0$xi.sigma2
    ctr$xi.omega <- 1/rgamma(1, 1, 1 + 1/ctr$omega)
    ctr$omega <- 1/rgamma(1, (ctr$sum.term + 1) / 2,
                          ctr$sum.f.exp / (2 * ctr$sigma2) + 1 / ctr$xi.omega)
    # ctr$omega <- 1/rgamma(1, (ctr$sum.term + 1) / 2,
    #                       ctr$sum.f.exp / (2 * ctr$sigma2) + 1 / 2)


    # if (any(is.nan(B0$gamma)) || is.nan(ctr$sigma2) || is.nan(ctr$xi.sigma2) ||
    #     is.nan(ctr$xi.omega) || is.nan(ctr$omega)) {
    #   cat(B0$gamma, ctr$sigma2, ctr$xi.sigma2, ctr$xi.omega, ctr$omega)
    #   stop("NaN values")
    # }

    # + Update tree split param ----
    # if (ctr$b > model$ctr$n.burn / 2 & model$ctr$split.prior == "logit") {
    #   alpha.new <- abs(rnorm(1, model$tree.params[1], model$tree.params[3]))
    #   alpha.MH <- do.call(sum, lapply(ctr$trees, function(tr) {
    #     tr$pTree.alpha(alpha.new, model$tree.params[1], 0)
    #   })) + dgamma(alpha.new, model$ctr$n.trees * model$tree.params[2],
    #                model$ctr$n.trees, log = TRUE) -
    #     dgamma(model$tree.params[1], model$ctr$n.trees * model$tree.params[2],
    #            model$ctr$n.trees, log = TRUE)
    #   if (lrunif[ctr$b] < alpha.MH) {
    #     model$tree.params[1] <- alpha.new
    #   }
    # }


    # + Record updates ----
    if (ctr$record) {
      if (model$intercept) {
        #intercept <- model$rescale$Y * B0$gamma[1] + attr(model$Y, "scaled:center") -
        #  sum((B0$gamma * attr(model$Z, "scaled:center"))[2:ctr$pZ])
        dgn$gamma[ctr$record,] <- B0$gamma #c(intercept, B0$gamma[2:ctr$pZ])
      } else {
        dgn$gamma[ctr$record,] <- B0$gamma
      }
      dgn$sigma2[ctr$record] <- ctr$sigma2
      dgn$fhat <- dgn$fhat + ctr$fhat
      if (model$ctr$diagnostics) {
        dgn$omega[ctr$record] <- ctr$omega
        dgn$tau[ctr$record,] <- ctr$tau
        # dgn$alpha[ctr$record] <- model$tree.params[1]
      }
    }

    # + Progress indicator ----
    if (model$ctr$verbose) {
      if (ctr$b%%100 == 0 & !model$ctr$par) {
        cat(".")
        # cat("nterm=", ctr$sum.term,
        #     ", mem/node", object_size(ctr$trees)/(ctr$sum.term * 1000^2),
        #     ", dgn = ", object_size(dgn), "\n")
        #gc(verbose = F);
      }
      if (ctr$b > model$ctr$n.burn &
          (ctr$b - model$ctr$n.burn) %% model$ctr$rec == 0) {
        if (model$ctr$par) cat("Core", model$ctr$core, ": ")
        cat((ctr$b - model$ctr$n.burn), "\n")
      }
      if (ctr$b == model$ctr$n.burn) {
        if (model$ctr$par) cat("Core", model$ctr$core, ": ")
        ctr$time <- proc.time()[3] - ctr$time
        if (ctr$time > 3600) {
          cat("\nBurn-in time:", round(ctr$time/3600, 2), "hours")
        } else if (ctr$time > 60) {
          cat("\nBurn-in time:", round(ctr$time/60, 2), "minutes")
        } else {
          cat("\nBurn-in time:", round(ctr$time, 2), "seconds")
        }
        ctr$time <- model$ctr$n.iter * ctr$time / model$ctr$n.burn
        if (ctr$time > 3600) {
          cat("\nEstimated time to completion:", round(ctr$time/3600, 2), "hours\n")
        } else if (ctr$time > 60) {
          cat("\nEstimated time to completion:", round(ctr$time/60, 2), "minutes\n")
        } else {
          cat("\nEstimated time to completion:", round(ctr$time, 2), "seconds\n")
        }
      }
    }
  }

  # ---- Prepare output ----
  DLMtemp <- lapply(1:length(dgn$DLM), function(i) do.call(rbind, dgn$DLM[[i]]))
  model$DLM <- as.data.frame(do.call(rbind, DLMtemp))
  colnames(model$DLM) <- c("Iter", "Tree", "xmin", "xmax", "tmin", "tmax", "est")
  model$DLM$est <- model$DLM$est * model$rescale$Y / model$Mo$Xscale
  if (length(model$Mo$Xsplits) == 0)
    model$Mo$X <- (model$Mo$X * model$Mo$Xscale) + model$Mo$Xmean
  rownames(model$DLM) <- NULL
  model$sigma2 <- dgn$sigma2 * model$rescale$Y^2
  if (model$intercept)
    dgn$gamma[,-1] <- sapply(2:ncol(dgn$gamma), function(i) dgn$gamma[,i] * model$rescale$Z[i])
  else
    dgn$gamma <- sapply(1:ncol(dgn$gamma), function(i) dgn$gamma[,i] * model$rescale$Z[i])
  model$gamma <- dgn$gamma
  colnames(model$gamma) <- colnames(model$Z)
  model$fhat <- model$rescale$Y * dgn$fhat /
    floor(model$ctr$n.iter/model$ctr$n.thin)
  if (model$ctr$diagnostics) {
    model$diagnostics <- list()
    model$diagnostics$omega <- dgn$omega
    model$diagnostics$tau <- dgn$tau
    # model$diagnostics$alpha <- dgn$alpha
    model$diagnostics$accept <- do.call(rbind.data.frame,
                                        lapply(1:length(dgn$accept),
                                               function(i) cbind.data.frame(dgn$accept[[i]])))
    model$diagnostics$term.nodes <- dgn$term.nodes
  }

  rm(dgn, ctr)
  invisible()
}
