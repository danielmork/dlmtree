predict.dlmtree <- function(object, new.data, new.exposure.data, ...,
                            ci.level = 0.95, type = "response", outcome = NULL,
                            fixed.idx = list(), est.dlm = FALSE, verbose = TRUE)
{
  `%notin%` <- Negate(`%in%`)
  out <- list()
  new.exposure.data <- matrix(new.exposure.data, ncol = object$pExp)
  if (!is.data.frame(new.data))
    step("`new.data` must be a data.frame with same colunm names as original model")
  if (!all(object$modNames %in% colnames(new.data)))
    stop("`new.data` must have the same colunm names as the original model")
  if (!is.numeric(new.exposure.data) || !is.matrix(new.exposure.data))
    stop("`new.exposure.data` must be a numeric matrix")
  if (nrow(new.data) != nrow(new.exposure.data))
    stop("`new.data` and `new.exposure.data` must have same number of rows")
  if (ncol(new.exposure.data) != object$pExp)
    stop("`new.exposure.data` must have the same number of measurments (columns) as original exposure data")

  ci.lims <- c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2)

  # ---- Setup covariate matrix and modifiers ----
  new.mf <- model.frame(delete.response(object$terms), new.data,
                        na.action = na.fail, drop.unused.levels = TRUE,
                        xlev = object$termLevels)
  z <-      model.matrix(delete.response(object$terms), new.mf,
                         xlev = object$termLevels)
  # z <-      z[, sort(object$QR$pivot[ seq_len(object$QR$rank) ]) ]
  # if (length(object$droppedCovar) > 0 & object$verbose)
  #   warning("variables {", paste0(object$droppedCovar, collapse = ", "),
  #           "} dropped from original due to perfect collinearity\n")
  n <-      nrow(z)
  mod <-    lapply(object$modNames, function(m) {
    if (!is.numeric(object$MoUnique[[m]])) {
      if (!all(unique(new.data[[m]]) %in% object$MoUnique[[m]]))
        stop("column ", m, " of `new.data` contains additional categories not in original dataset")
    } else {
      if (min(new.data[[m]]) < object$MoUnique[[m]][1] ||
          max(new.data[[m]]) > object$MoUnique[[m]][2])
        if (verbose)
          warning("column ", m, "of `new.data` exceeds range of original dataset")
    }
    new.data[[m]]
  })
  names(mod) <- object$modNames

  # ---- Predict z^T * gamma ----
  ztg.draws <- z %*% t(object$gamma)
  out$ztg <- rowMeans(ztg.draws)
  out$ztg.lims <- apply(ztg.draws, 1, quantile, probs = ci.lims)

  # ---- Predict DLMs ----
  mark <- ceiling(nrow(object$TreeStructs) / 42)
  if (verbose)
    cat(paste0("Reanalyzing trees for new.data: % complete\n",
               "[0--------25--------50--------75--------100]\n '"))

  draws <- lapply(1:object$mcmcIter, function(i) matrix(0.0, n, object$pExp))

  if (is.null(object$fixedIdx)) {
    for (i in 1:nrow(object$TreeStructs)) {
      if (verbose && (i %% mark == 0))
        cat("'")
      Iter <- object$TreeStructs$Iter[i]
      rule <- object$TreeStructs$Rule[i]
      if (rule == "")
        idx <- 1:n
      else
        idx <- which(eval(parse(text = rule)))

      if (object$dlmType == "gp") {
        est <- matrix(rep(as.numeric(object$TreeStructs[i, -c(1:4)]), each = length(idx)), length(idx), object$pExp)
        draws[[Iter]][idx,] <- draws[[Iter]][idx,] + est
      } else {
        t <- object$TreeStructs$tmin[i]:object$TreeStructs$tmax[i]
        est <- object$TreeStructs$est[i]
        draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
      }
    }
  } else {
    for (i in 1:nrow(object$TreeStructs)) {
      if (i %% mark == 0)
        cat("'")
      Iter <- object$TreeStructs$Iter[i]
      idx <- fixed.idx[[object$TreeStructs$fixedIdx[i] + 1]]

      if (object$dlmType == "gp") {
        est <- matrix(rep(as.numeric(object$TreeStructs[i, -c(1:4)]), each = length(idx)), length(idx), object$pExp)
        draws[[Iter]][idx,] <- draws[[Iter]][idx,] + est
      } else {
        t <- object$TreeStructs$tmin[i]:object$TreeStructs$tmax[i]
        est <- object$TreeStructs$est[i]
        draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
      }
    }
  }


  draws <- array(do.call(c, draws), c(n, object$pExp, object$mcmcIter))

  if (est.dlm) {
    if (verbose)
      cat("\nestimating individualized DLMs...")
    out$dlmest <- sapply(1:object$pExp, function(t) {
      rowMeans(draws[,t,,drop=F])
    })
    out$dlmest.lower <- sapply(1:object$pExp, function(t) {
      apply(draws[,t,,drop=F], 1, quantile, probs = 0.025)
    })
    out$dlmest.upper <- sapply(1:object$pExp, function(t) {
      apply(draws[,t,,drop=F], 1, quantile, probs = 0.975)
    })
  }
  if (verbose)
    cat("\nestimating exposure effects...")
  fhat.draws <- do.call(rbind, lapply(1:n, function(i) {
    t(t(draws[i,,]) %*% new.exposure.data[i,])
  }))
  out$fhat <- rowMeans(fhat.draws)
  out$fhat.lims <- apply(fhat.draws, 1, quantile, probs = ci.lims)

  # ---- Outcome predictions ----
  if (type == "response") {
    y.draws <- ztg.draws + fhat.draws +
      sapply(1:object$mcmcIter, function(i) rnorm(nrow(ztg.draws), 0, sqrt(object$sigma2[i])))
    out$y <- rowMeans(y.draws)
    out$y.lims <- apply(y.draws, 1, quantile, probs = ci.lims)

    return(out)
  } else if (type == "waic") {
    if (length(outcome) != nrow(ztg.draws))
      stop("must provide compplete data outcome to calculate WAIC")
    probs <- sapply(1:object$mcmcIter, function(i) dnorm(outcome, ztg.draws[,i] + fhat.draws[,i], sqrt(object$sigma2[i])))
    LPD <- sum(log(rowMeans(probs)))
    pwaic <- sum(apply(log(probs), 1, var))
    return(list(waic = -2 * (LPD - pwaic), LPD = LPD, pwaic = pwaic))
  }
}
