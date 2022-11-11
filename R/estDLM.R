estDLM <- function(object,
                   new.data,
                   group.index,
                   conf.level = 0.95,
                   cenval = 0,
                   return.mcmc = FALSE,
                   mem.safe = FALSE,
                   verbose = TRUE)
{
  if (class(object) != 'dlmtree')
    stop("`object` must be return from model `dlmtree`")
  if (!all(object$modNames %in% colnames(new.data)))
    stop("`new.data` must have the same colunm names as the original model")

  out <- list()
  class(out) <- "estDLM"
  out$conf.level <- conf.level
  ci.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  # Load modifier data
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

  # Name group indices
  if (!is.list(group.index))
    group.index <- list(group.index)
  if (is.null(names(group.index)))
    names(group.index) <- 1:length(group.index)
  if (any(sapply(group.index, function(x) any(x > nrow(new.data)))))
    stop("at least one item of `group.index` exceeds range of data")
  out$n <- lapply(group.index, length)
  out$groupIndex <- group.index

  tempMod <- data.frame(Rule = unique(object$TreeStructs$Rule))
  rules <- strsplit(tempMod$Rule, " & ", T)
  mark <- ceiling(length(rules) / 42)

  # Analyze trees for each group index list
  for (i in 1:length(group.index)) {
    if (names(group.index)[i] == "")
      names(group.index[i]) <- i
    if (verbose)
      cat(paste0("Reanalyzing trees for group ", names(group.index)[i], "\n",
                 "% complete\n[0--------25--------50--------75--------100]\n '"))
    modDat <- lapply(mod, function(k) k[group.index[[i]]])
    names(modDat) <- names(mod)
    ri <- ruleIdx(modDat, mem.safe)
    tempMod[[paste0("Weight", i)]] <- sapply(1:length(rules), function(i) {
      if (verbose && (i %% mark == 0))
        cat("'")
      ri$returnWeight(rules[[i]])
    })
    rm(ri)
    if (verbose)
      cat("\n")
  }

  if (verbose)
    cat("Calcuating DLMs...")
  DLM <- merge.data.frame(object$TreeStructs, tempMod, by = "Rule")
  if (return.mcmc)
    out$mcmc <- list()
  out$dlmMean <- list()
  out$dlmCI <- list()
  out$dlmCum <- list()


  # ---- dlmType: GP ----
  if (object$dlmType == "gp") {
    for (i in 1:length(group.index)) {
      gDLM <- DLM[,5:(4 + object$pExp)] * DLM[[paste0("Weight", i)]]
      mcmc <-
        sapply(1:max(DLM$Iter), function(i) {
          rowSums(
            sapply(1:max(DLM$Tree), function(t) {
              colSums(gDLM[which(DLM$Iter == i & DLM$Tree == t),, drop = F])
            })
          )
        })
      if (return.mcmc)
        out$mcmc[[names(group.index)[i]]] <- mcmc
      out$dlmMean[[names(group.index)[i]]] <- rowMeans(mcmc)
      out$dlmCI[[names(group.index)[i]]] <- apply(mcmc, 1, quantile, probs = ci.lims)
      out$dlmCum[[names(group.index)[i]]] <- c(mean = mean(rowMeans(mcmc)), quantile(colSums(mcmc), probs = ci.lims))
    }
    out$dlFunction <- "dlm"


  # ---- dlmType: TDLM ----
  } else if (object$dlFunction == "dlm") {
    for (i in 1:length(group.index)) {
      DLM$w.est <- DLM$est * DLM[[paste0("Weight", i)]]
      mcmc <- dlmEst(as.matrix(DLM[,c("Iter", "Tree", "tmin", "tmax", "w.est")]),
                     object$pExp, object$mcmcIter)
      if (return.mcmc)
        out$mcmc[[names(group.index)[i]]] <- mcmc
      out$dlmMean[[names(group.index)[i]]] <- rowMeans(mcmc)
      out$dlmCI[[names(group.index)[i]]] <- apply(mcmc, 1, quantile, probs = ci.lims)
      out$dlmCum[[names(group.index)[i]]] <- c(mean = mean(rowMeans(mcmc)), quantile(colSums(mcmc), probs = ci.lims))
    }
    out$dlFunction <- "dlm"

  # ---- dlmType: TDLNM ----
  } #else {
  #   for (i in 1:length(group.index)) {
  #     DLM$w.est <- DLM$est * DLM[[paste0("Weight", i)]]
  #     if (is.na(object$SE[1])) {
  #       cen.quant <- which.min(abs(object$Xsplits - cenval))
  #       out$mcmc[[names(group.index)[i]]] <-
  #         dlnmEst(as.matrix(DLM[,c("Iter", "Tree", "xmin", "xmax",
  #                                  "tmin", "tmax", "w.est")]),
  #                 object$Xsplits, object$pExp,  object$mcmcIter, cen.quant, 0)
  #     } else {
  #       out$mcmc[[names(group.index)[i]]] <-
  #         dlnmEst(as.matrix(DLM[,c("Iter", "Tree", "xmin", "xmax",
  #                                  "tmin", "tmax", "w.est")]),
  #                 object$Xsplits, object$pExp,  object$mcmcIter,
  #                 cenval, mean(object$SE))
  #     }
  #   }
  #   out$dlFunction <- "dlnm"
  # }
  lags <- length(out$dlmMean[[1]])
  out$plotData <- do.call(rbind, lapply(names(group.index), function(n) {
    data.frame(group = n, time = 1:lags, est = out$dlmMean[[n]], 
               lower = out$dlmCI[[n]][1,], upper = out$dlmCI[[n]][2,])
  }))

  return(out)
}



ruleIdx <- function(mod, mem.safe = FALSE)
{
  self <- environment()
  n <- length(mod[[1]])
  idxList <- list()
  `%notin%` <- Negate(`%in%`)

  returnWeight = function(rules)
  {
    if (length(rules) == 0)
      return(1)

    if (mem.safe) {
      return(length(Reduce(cppIntersection, lapply(rules, function(i) which(eval(parse(text = i))))))/self$n)
    } else {
      missingIdx <- sapply(rules, function(i) is.null(idxList[[i]]))
      newIdx <- list()
      idx <- list()
      if (any(missingIdx)) {
        newIdx <- lapply(rules[missingIdx], function(i) which(eval(parse(text = i))))
        names(newIdx) <- rules[missingIdx]
        idxList <<- append(idxList, newIdx)
      }
      return(length(Reduce(cppIntersection, lapply(rules, function(i) idxList[[i]])))/self$n)
    }
  }

  return(self)
}
