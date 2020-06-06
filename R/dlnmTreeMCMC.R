#' dlnmTreeMCMC
#'
#' @param model Environment with model data and control parameters
#' @param ctr Environment with MCMC control parameters
#' @param dgn Environment for storing output and diagnostics
#'
#' @description Updates all trees in model, alters environments, no return.
#'
#' @return invisible()
#'
dlnmTreeMCMC <- function(model, ctr, dgn)
{
  n.trees <- model$ctr$n.trees
  ctr$R <- ctr$R + ctr$Rmat[[1]]
  ctr$sum.term <- 0
  ctr$sum.f.exp <- 0
  DLM <- list()

  # Sample tree update steps
  step <- sample.int(4, n.trees, replace = T, prob = model$step.prob)
  lrunif <- log(runif(n.trees))

  # Create new tree proposals
  new.trees <- lapply(1:model$ctr$n.trees, function(i) ctr$trees[[i]]$copy())
  update.proposals <- lapply(1:model$ctr$n.trees, function(i) {
    dlnmTreeUpdate(new.trees[[i]], step[i], model)
  })

  # MCMC on trees
  for (t in 1:n.trees) {
    proposal <- update.proposals[[t]]
    tree.var <- ctr$tau[t] * ctr$omega
    new.trees[[t]]$update <- TRUE

    # Combine terminal node counts
    old.term <- ctr$trees[[t]]$list.terminal()
    Xd.old <- do.call(cbind, lapply(old.term, function(i) i$get.counts(model, ctr)))
    new.term <- new.trees[[t]]$list.terminal()
    Xd.new <- do.call(cbind, lapply(new.term, function(i) i$get.counts(model, ctr)))

    # Gather data concerning old tree
    if (length(old.term) == 1) {
      orig <- dlnmMHR1(ctr$R, model$Z, ctr$XX1, ctr$ZtX1, ctr$VgZtX1,
                       ctr$Vg, ctr$X1, c(0), 1 / tree.var, sqrt(ctr$sigma2))
    } else {
      ZtX.old <- do.call(cbind, lapply(old.term, function(i) i$ZtX))
      VgZtX.old <- do.call(cbind, lapply(old.term, function(i) i$VgZtX))
      orig <- dlnmMHR(ctr$R, model$Z, ZtX.old, VgZtX.old,
                      ctr$Vg, Xd.old, c(0), 1 / tree.var, sqrt(ctr$sigma2))
    }

    if (is.nan(orig$VThetaLogDet)) {
      cat("Iter: ", ctr$b, "Tree: ", t, "\n")
      cat("omega:", ctr$omega, "tau_a:", ctr$tau[t], "\n")
      stop("NaN MHR")
    }

    # Propose new tree
    if (proposal$Accept & ncol(Xd.new) == length(new.term)) {
      if (length(new.term) == 1) {
        new <- dlnmMHR1(ctr$R, model$Z, ctr$XX1, ctr$ZtX1, ctr$VgZtX1,
                        ctr$Vg, ctr$X1, orig$ZtY, 1 / tree.var, sqrt(ctr$sigma2))
      } else {
        ZtX.new <- do.call(cbind, lapply(new.term, function(i) i$ZtX))
        VgZtX.new <- do.call(cbind, lapply(new.term, function(i) i$VgZtX))
        new <- dlnmMHR(ctr$R, model$Z, ZtX.new, VgZtX.new,
                       ctr$Vg, Xd.new, orig$ZtY,
                       1 / tree.var, sqrt(ctr$sigma2))
      }

      tree.mhr <- proposal$MHR + (new$VThetaLogDet - orig$VThetaLogDet) -
        (ctr$n + 1)/2 * (log(new$Beta/2 + 1/ctr$xi.sigma2) - log(orig$Beta/2 + 1/ctr$xi.sigma2)) -
        log(tree.var) * (length(new.term) - length(old.term))/2

      # proposal$MHRV <- (new$VThetaLogDet - orig$VThetaLogDet)
      # proposal$MHRB <- -(ctr$n + 1)/2 * (log(new$Beta/2 + 1/ctr$xi.sigma2) - log(orig$Beta/2 + 1/ctr$xi.sigma2))
      # proposal$MHRD <- -log(tree.var) * (length(new.term) - length(old.term))/2
      # proposal$tree <- t

      if (is.nan(tree.mhr)) {
        tree.mhr <- -Inf
        cat("NaN produced in MHR, update skipped. Details:")
        cat("Iter: ", ctr$b, "Tree: ", t, "NB:", new$Beta, "OB:", orig$Beta, "\n")
        cat("omega:", ctr$omega, "tau_a:", ctr$tau[t], "\n")
      } else {
        proposal$MHR2 <- tree.mhr
      }
    } else {
      # proposal$MHRV <- 0
      # proposal$MHRB <- 0
      # proposal$MHRD <- 0
      # proposal$tree <- t
      tree.mhr <- -Inf
    }

    # Accept/reject new tree
    if (lrunif[t] < tree.mhr) {
      proposal$Accept <- 2
      ctr$trees[[t]] <- new.trees[[t]]$copy()
      term <- new.term
      draw <- new$ThetaDraw
      ctr$TLiT[t] <- new$TLiT
      ctr$Rmat[[t]] <- new$Yhat
      ctr$n.term[t] <- length(new.term)
    } else {
      term <- old.term
      draw <- orig$ThetaDraw
      ctr$TLiT[t] <- orig$TLiT
      ctr$Rmat[[t]] <- orig$Yhat
    }

    # Update with new draw
    ctr$xi.tau[t] <- 1/rgamma(1, 1, 1 + 1 / (ctr$tau[t]))
    ctr$tau[t] <- 1/rgamma(1, (length(term) + 1) / 2,
                           (ctr$TLiT[t] / (2 * ctr$sigma2 * ctr$omega)) +
                             (1 / ctr$xi.tau[t]))
    # ctr$tau[t] <- 1/rgamma(1, (length(term) + 1) / 2,
    #                        (ctr$TLiT[t] / (2 * ctr$sigma2 * ctr$omega)) +
    #                          (1 / 2))
    ctr$sum.term <- ctr$sum.term + length(term)
    ctr$sum.f.exp <- ctr$sum.f.exp + ctr$TLiT[t] / ctr$tau[t]

    if (t < n.trees) {
      ctr$R <- ctr$R - ctr$Rmat[[t]] + ctr$Rmat[[t+1]]
    }

    # Save every kth draw after burn-in
    if (ctr$record) {
      rules <- do.call(rbind, lapply(term, function(i) c(i$xmin, i$xmax, i$tmin, i$tmax)))
      # dgn$DLM[[(ctr$record - 1) * n.trees + t]] <- cbind(ctr$record, t, rules, draw)
      DLM[[t]] <- cbind(ctr$record, t, rules, draw)
      if (model$ctr$diagnostics) {
        proposal$Var <- ctr$sigma2 * tree.var
        dgn$accept[[(ctr$record - 1) * n.trees + t]] <- proposal
      }
    }

  }
  if (ctr$record) {
    dgn$DLM[[ctr$record]] <- DLM
    if (model$ctr$diagnostics) {
      dgn$term.nodes[ctr$record, ] <- ctr$n.term
    }
  }
  #cat(ctr$b, ctr$sigma2, ctr$omega, paste0(which.max(ctr$tau), ": ", max(ctr$tau)),"\n")
  invisible()
}
