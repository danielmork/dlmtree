dlmtree <- function(formula,
                    data,
                    exposure.data,
                    ...,
                    tree.modifiers = "all",
                    modifier.splits = 20,
                    dlm.type = "nested",
                    #ver = 1,
                    covariance.type = "exponential",
                    # exposure.splits = 0,
                    exposure.center = "median",
                    # exposure.se = NULL,
                    n.trees = 20,
                    n.burn = 2000,
                    n.iter = 5000,
                    n.thin = 10,
                    family = "gaussian",
                    binomial.size = 1,
                    fixed.tree.idx = NULL,
                    tree.params.mod = c(.95, 2),
                    tree.params.tdlm = c(.95, 2),
                    step.prob.mod = c(.25, .4, .1),
                    step.prob.tdlm = c(.25, .25),
                    zeta.prior = 0.5,
                    shrinkage = 1,
                    subset = 1:nrow(data),
                    save.data = TRUE,
                    verbose = TRUE,
                    diagnostics = FALSE,
                    initial.params = NULL)
{
  options(stringsAsFactors = FALSE)
  model <- list()

  # ---- Check inputs ----
  if (!is.data.frame(data))
    stop("`data` must be a data.frame")
  if (!is.numeric(exposure.data))
    stop("`exposure.data` must be a numeric matrix")
  if (nrow(data) != nrow(exposure.data))
    stop("`data` and `exposure.data` must have same number of observations")
  # if (!is.null(exposure.se)) {
  #   if (!is.numeric(exposure.se))
  #     stop("`exposure.se` must be a scalar or numeric matrix")
  #   if (length(exposure.se) == 1)
  #     exposure.se <- matrix(exposure.se, nrow(exposure.data), ncol(exposure.data))
  #   if (!all(dim(exposure.data) == dim(exposure.se)))
  #     stop("`exposure.se` and `exposure.data` must have same dimensions")
  # }

  # iteration control
  if (n.iter < n.thin * 10)
    stop("after thinning you will be left with less than 10 MCMC samples,",
         " increase the number of iterations!")
  if (!(dlm.type %in% c("gp", "shared", "nested", "tdlm", "tdlm2")))
    stop("`dlm.type` must be one of `gp`, `shared`, or `nested`")

  # response type
  if (!(family %in% c("gaussian", "logit")))
    stop("`family` must be one of `gaussian`, or `logit`")
  model$binomial <- FALSE
  if (family == "logit") {
    model$binomial <- TRUE
    if (length(binomial.size) == 1)
      binomial.size <- rep(binomial.size, nrow(data))
    if (length(binomial.size) != nrow(data))
      stop("`binomial.size` must be positive integer and same length as data")
    model$binomialSize <- force(binomial.size)
  }

  # tree parameters
  if (length(tree.params.mod) != 2 || length(tree.params.tdlm) != 2) {
    stop("tree.params.* must have length 2")
  } else {
    if (tree.params.mod[1] > 1 || tree.params.mod[1] < 0 ||
        tree.params.tdlm[1] > 1 || tree.params.tdlm[1] < 0)
      stop("tree.params.*[1] must be between 0 and 1")
    if (tree.params.mod[2] < 0 || tree.params.tdlm[2] < 0)
      stop("tree.params.*[2] must be greater than 0")
  }

  # step probabilities
  if (any(step.prob.mod < 0) || any(step.prob.mod > 1) ||
      any(step.prob.tdlm < 0) || any(step.prob.tdlm > 1))
    stop("step.prob.* components must be between zero and one")


  # ---- Model control arguments ----
  model$nTrees <-   n.trees
  model$nBurn <-    n.burn
  model$nIter <-    n.iter
  model$nThin <-    n.thin
  model$mcmcIter <- floor(n.iter / n.thin)
  model$dlmType <-  dlm.type
  model$family <-   family
  model$zeta <-     zeta.prior
  model$shrinkage <-     shrinkage
  model$verbose <-       verbose
  model$diagnostics <-   diagnostics
  model$treePriorMod <-  tree.params.mod
  model$treePriorTDLM <- tree.params.tdlm
  model$stepProbMod <-   c(step.prob.mod[1], step.prob.mod[1],
                           step.prob.mod[2], step.prob.mod[3])
  model$stepProbMod <-   model$stepProbMod / sum(model$stepProbMod)
  model$stepProbTDLM <-  c(step.prob.tdlm[1], step.prob.tdlm[1],
                           step.prob.tdlm[2])
  model$stepProbTDLM <-  model$stepProbTDLM / sum(model$stepProbTDLM)

  if (model$verbose)
    cat("Preparing data...\n")

  # ---- Create data subset ----
  # if (is.null(subset))
  #   subset <- 1:nrow(data)
  # else
  if (!is.integer(subset) || any(subset < 1) || any(subset > nrow(data)))
    stop("invalid subset, must be integers within range of data length")
  data <- data[subset,]
  exposure.data <- exposure.data[subset,]

  # ---- Setup fixed effect model ----
  model$formula <-    as.formula(formula)
  model$terms <-      terms.formula(model$formula, data = data)
  if (!attr(model$terms, "response"))
    stop("no valid response in formula")
  if (length(which(attr(model$terms, "term.labels") %in% colnames(data))) == 0)
    stop("no valid variables in formula")
  mf <-               model.frame(model$terms, data = data,
                                  na.action = NULL,
                                  drop.unused.levels = TRUE)
  model$termLevels <- .getXlevels(delete.response(model$terms), mf)
  model$intercept <-  ifelse(attr(model$terms, "intercept"), TRUE, FALSE)
  if (any(is.na(mf)))
    stop("missing values in model data, use `complete.cases()` to subset data")

  # ---- Drop collinear variables ----
  model$Y <-   model.response(mf)
  model$Z <-   model.matrix(model$terms, data = mf)
  # model$QR <-        qr(crossprod(model$Z))
  model$Znames <- colnames(model$Z)#[sort(QR$pivot[seq_len(QR$rank)])]
  model$droppedCovar <- c()#colnames(model$Z)[ model$QR$pivot[ -seq_len(model$QR$rank) ] ]
  # model$Z <-   model$Z[, sort(model$QR$pivot[ seq_len(model$QR$rank) ]) ]
  # model$Z <-   scaleModelMatrix(model$Z)
  # if (length(model$droppedCovar) > 0 & model$verbose)
  #   warning("variables {", paste0(model$droppedCovar, collapse = ", "),
  #           "} dropped due to perfect collinearity\n")

  # ---- Scale data ----
  if (model$family == "gaussian") {
    model$Ymean <-  sum(range(model$Y))/2
    model$Yscale <- diff(range(model$Y - model$Ymean))
    model$Y <-      (model$Y - model$Ymean) / model$Yscale
  } else {
    model$Yscale <- 1
    model$Ymean <-  0
    model$Y <-      scale(model$Y, center = 0, scale = 1)
  }
  model$Y <-       c(model$Y)
  model$Zscale <-  rep(1, ncol(model$Z))#attr(model$Z, "scaled:scale")
  model$Zmean <-   rep(0, ncol(model$Z))#attr(model$Z, "scaled:center")
  model$Z <-       matrix(model$Z, nrow(model$Z), ncol(model$Z))


  model$initParams <- rep(0, ncol(model$Z))
  if (!is.null(initial.params)) {
    names(model$initParams) <- colnames(model$Z)
    if (sum(names(initial.params) %in% colnames(model$Z)) > 0) {
      na <- names(initial.params[ # get matching names 
        which(names(initial.params) %in% colnames(model$Z))])
      model$initParams[na] <- initial.params[na]
    }
  }


  # ---- Exposure splits ----
  model$pExp <-     ncol(exposure.data)
  model$timeProb <- rep(1 / (model$pExp - 1), model$pExp - 1)
  model$X <-        exposure.data
  model$Xrange <-   range(model$X)
  # if (is.null(exposure.se)) {
    model$smooth <- FALSE
    model$SE <-     matrix(0.0, 0, 0)
  # } else {
  #   model$smooth <- TRUE
  #   model$SE <-     exposure.se[subset, ]
  # }

  # if (length(exposure.splits) == 1) {
    # no splits -- DLM
    # if (exposure.splits == 0) {
      model$dlFunction <- "dlm"
      model$splitProb <-  as.double(c())
      model$Xsplits <-    as.double(c())
      model$nSplits <-    0
      model$Xscale <-     sd(model$X)
      model$X <-          model$X / model$Xscale
      model$Tcalc <-      sapply(1:ncol(model$X), function(i) rowSums(model$X[, 1:i, drop = F]))

      # splits defined by quantiles of exposure
    # } else {
    #   model$dlFunction <- "dlnm"
    #   if (is.list(exposure.splits)) {
    #     stop("exposure.splits must be a scalar or list with two inputs: 'type' and 'split.vals'")
    #   } else {
    #     model$Xsplits <- force(sort(unique(quantile(model$X,
    #                                                 (1:(exposure.splits - 1)) /
    #                                                   exposure.splits))))
    #     model$nSplits <- force(length(model$Xsplits))
    #     model$splitProb <- force(rep(1 / model$nSplits, model$nSplits))
    #   }
    # }

  # splits defined by specific values or quantiles
  # } else {

  #   model$dlFunction <- "dlnm"
  #
  #   # if exposure.splits entered incorrectly, infer user input and inform
  #   if (!is.list(exposure.splits)) {
  #     if (any(exposure.splits > 1 | exposure.splits < 0)) {
  #       cat("exposure.splits entered as numeric vector, assuming values are exposure splitting points\n")
  #       exposure.splits <- list("type" = "values",
  #                               "split.vals" = exposure.splits)
  #     } else {
  #       cat("exposure.splits entered as numeric vector, assuming values are exposure splitting quantiles\n")
  #       exposure.splits <- list("type" = "quantiles", "split.vals" = exposure.splits)
  #     }
  #   }
  #
  #   # use specific values as splitting points
  #   if (exposure.splits$type == "values") {
  #     model$Xsplits <- force(sort(unique(exposure.splits$split.vals)))
  #     model$Xsplits <- force(model$Xsplits[which(model$Xsplits > min(model$X) &
  #                                                  model$Xsplits < max(model$X))])
  #
  #     # use specific quantiles as splitting points
  #   } else if (exposure.splits$type == "quantiles") {
  #     if (any(exposure.splits$split.vals > 1 | exposure.splits$split.vals < 0))
  #       stop("`exposure.splits$split.vals` must be between zero and one if using quantiles")
  #     model$Xsplits <- force(sort(unique(quantile(model$X,
  #                                                 exposure.splits$split.vals))))
  #     model$Xsplits <- force(model$Xsplits[which(model$Xsplits > min(model$X) &
  #                                                  model$Xsplits < max(model$X))])
  #   } else {
  #     stop("`exposure.splits$type` must be one of `values` or `quantiles`")
  #   }
  #
  #   model$nSplits <- force(length(model$Xsplits))
  #   if (model$nSplits == 0)
  #     stop("no exposure splits specified, please check `exposure.splits` input")
  #
  #   model$splitProb <- force(rep(1 / model$nSplits, model$nSplits))
  # }


  # Precalculate counts below each splitting values
  # if (length(model$Xsplits) > 0) {
  #   model$Xscale <- 1
  #   model$Tcalc <- force(sapply(1:ncol(model$X), function(i) {
  #     rep(i / model$Xscale, nrow(model$X)) }))
  #   if (model$smooth) {
  #     model$Xcalc <- force(sapply(model$Xsplits, function(i) {
  #       rowSums(pnorm((i - model$X) / model$SE)) }) / model$Xscale)
  #   } else {
  #     model$Xcalc <- force(sapply(model$Xsplits, function(i) {
  #       rowSums(model$X < i) }) / model$Xscale)
  #   }
  #   if (exposure.center == "median")
  #     model$Xcenter <- median(model$X)
  #   else
  #     model$Xcenter <- exposure.center
  #   model$XcenterIdx <- mean(c(min(which(model$Xsplits - model$Xcenter > 0)),
  #                              max(which(model$Xsplits - model$Xcenter <= 0))))
  # }


  # ---- Setup modifier variables ----
  model$modNames <- tree.modifiers
  if (length(model$modNames) == 1)
    if (model$modNames == "all")
      model$modNames <- colnames(data)[-which(colnames(data) == all.vars(model$formula[[2]]))]

  model$Mo <- lapply(model$modNames, function(m) {
    if (!(m %in% colnames(data)))
      stop("one or more `tree.modifiers` specified is not a column of `data`")
    data[[m]]
  })
  model$MoUnique <- lapply(model$modNames, function(m) {
    if (!is.numeric(data[[m]]))
      sort(unique(data[[m]]))
    else
      range(data[[m]])
  })
  names(model$Mo) <- names(model$MoUnique) <- model$modNames
  model$pM <- length(model$Mo)

  # Use quantiles of modifiers as possible splits
  model$modIsNum <- sapply(model$Mo, function(i) (is.numeric(i) || is.logical(i)))
  model$modSplitValRef <- list()
  model$modSplitValIdx <- list()
  model$modSplitIdx <- list()
  model$fullIdx <- 0:(nrow(data) - 1)
  for (i in 1:model$pM) {

    if (model$modIsNum[i]) {

      if (length(unique(model$Mo[[i]])) < modifier.splits) {
        uniqueVals <- sort(unique(model$Mo[[i]]))
        model$modSplitValRef[[i]] <-
          rowMeans(cbind(uniqueVals[-length(uniqueVals)], uniqueVals[-1]))

      } else {
        uniqueVals <- sort(unique(quantile(
          model$Mo[[i]], 1:modifier.splits/(modifier.splits + 1))))
        model$modSplitValRef[[i]] <-
          uniqueVals[which(uniqueVals > min(model$Mo[[i]]) &
                             uniqueVals < max(model$Mo[[i]]))]
      }
      model$modSplitValIdx[[i]] <- 0:(length(model$modSplitValRef[[i]]) - 1)
      model$modSplitIdx[[i]] <-
        lapply(model$modSplitValRef[[i]], function(j) which(model$Mo[[i]] < j) - 1)

    } else {
      if (length(unique(model$Mo[[i]])) != length(unique(model$Mo[[i]])))
        warning("one or more modifier categories not present in data subset")
      model$modSplitValRef[[i]] <- sort(unique(model$Mo[[i]]))
      model$modSplitValIdx[[i]] <- 0:(length(model$modSplitValRef[[i]]) - 1)
      model$modSplitIdx[[i]] <-
        lapply(model$modSplitValRef[[i]], function(j) which(model$Mo[[i]] == j) - 1)
    }
  }


  # ---- Run model ----
  if (is.null(fixed.tree.idx)) {
    if (model$dlmType %in% c("tdlm", "shared"))
      out <- dlmtreeTDLMGaussian(model)
    else if (model$dlmType %in% c("tdlm2", "nested"))
      # if (ver == 2)
        out <- dlmtreeTDLM_cpp(model)
      # else if (ver == 1)
      #   out <- dlmtreeTDLMNestedGaussian(model)
    else if (model$dlmType == "gp") {
      model$covarianceType = 0
      if (covariance.type == "exponential") {
        model$covarianceType = 1
        model$DistMat <- -toeplitz(0:(model$pExp - 1))
      } else if (covariance.type == "gaussian") {
        model$covarianceType = 1
        model$DistMat <- -toeplitz((0:(model$pExp - 1))^2)
      } else if (covariance.type == "identity") {
        model$covarianceType = 0
        model$DistMat <- diag(model$pExp)
      }
      out <- dlmtreeGPGaussian(model)
    }

  } else {
    idx <- sort(unique(do.call(c, fixed.tree.idx)))
    if (length(idx) != length(model$Y) | any(idx > length(model$Y)))
      stop("fixed.tree.idx must contain all indices once")
    model$fixedIdx <- lapply(fixed.tree.idx, function(i) i - 1)

    if (model$dlmType == "gp") {
      model$DistMat <- -toeplitz(0:(model$pExp - 1))
      out <- dlmtreeGPFixedGaussian(model)
    } else {
      out <- dlmtreeTDLMFixedGaussian(model)
    }
  }


  if (verbose)
    cat("\nCompiling results...\n")

  for (n in names(out))
    model[[n]] <- out[[n]]


  # ---- Prepare output ----
  model$Y <- model$Y * model$Yscale + model$Ymean
  model$fhat <- model$fhat * model$Yscale
  model$sigma2 <- model$sigma2 * (model$Yscale ^ 2)

  # rescale fixed effect estimates
  model$gamma <- sapply(1:ncol(model$gamma), function(i) {
    model$gamma[,i] * model$Yscale / model$Zscale[i] })
  if (model$intercept) {
    model$gamma[,1] <- model$gamma[,1] + model$Ymean
    if (ncol(model$Z) > 1)
      model$gamma[,1] <- model$gamma[,1] - model$gamma[,-1,drop=FALSE] %*% model$Zmean[-1]
  }
  colnames(model$gamma) <- model$Znames

  if (is.null(fixed.tree.idx))
    colnames(model$modProb) <- colnames(model$modCount) <-
      colnames(model$modInf) <- names(model$Mo)

  if (is.null(fixed.tree.idx)) {
    modNames <- names(model$Mo)
    splitRules <- strsplit(model$termRules, "&", T)
    rule <- sapply(splitRules, function(str) {
      paste0(lapply(sort(str), function(rule) {
        # no rule
        if (length(rule) == 0) {
          return("")
        # >=
        } else if (length(spl <- strsplit(rule, ">=", T)[[1]]) == 2) {
          return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1], "']] >= ",
                        model$modSplitValRef[[as.numeric(spl[1]) + 1]][
                          as.numeric(spl[2]) + 1]))
        # >
        } else if (length(spl <- strsplit(rule, "<", T)[[1]]) == 2) {
          return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1], "']] < ",
                        model$modSplitValRef[[as.numeric(spl[1]) + 1]][
                          as.numeric(spl[2]) + 1]))
        # in
        } else if (length(spl <- strsplit(rule, "[]", T)[[1]]) == 2) {
          inList <- paste0("c('", paste0(model$modSplitValRef[[
            as.numeric(spl[1]) + 1]][
              eval(parse(text = paste0("c(", spl[2], ")"))) + 1
            ], collapse = "','"), "')")
          return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1],
                        "']] %in% ", inList))
        # not in
        } else if (length(spl <- strsplit(rule, "][", T)[[1]]) == 2) {
          inList <- paste0("c('", paste0(model$modSplitValRef[[
            as.numeric(spl[1]) + 1]][
              eval(parse(text = paste0("c(", spl[2], ")"))) + 1
            ], collapse = "','"), "')")
          return(paste0("mod[['", modNames[as.numeric(spl[1]) + 1],
                        "']] %notin% ", inList))
        } else {
          return("")
        }
      }), collapse = " & ")})
    model$modPairs <- sort(table(do.call(c, lapply(splitRules, function(r) {
      if (length(r) == 0)
        return(NA)
      m <- sapply(strsplit(r, ">=|<|\\[\\]|\\]\\["), function(i) as.numeric(i[1]))
      if (length(m) < 2)
        return(modNames[m + 1])
      c <- combn(length(m), 2)
      return(unique(sapply(1:ncol(c), function(i) paste0(modNames[sort(m[c[,i]]) + 1], collapse = "-"))))
    })))) / (model$nTrees * model$mcmcIter)

    model$TreeStructs <- cbind.data.frame(rule, model$TreeStructs)
    model$modSplitIdx <- NULL
    model$fullIdx <- NULL

    if (model$dlmType == "gp") {
      colnames(model$TreeStructs) <- c("Rule", "Iter", "Tree", "modTerm",
                                       paste0("Lag", 1:model$pExp))
      model$TreeStructs[,-c(1:4)] <- model$TreeStructs[,-c(1:4)] *
        model$Yscale / model$Xscale
    } else {
      colnames(model$TreeStructs) <- c("Rule", "Iter", "Tree", "modTerm", "dlnmTerm",
                                       "xmin", "xmax", "tmin", "tmax", "est")
      model$TreeStructs$est <- model$TreeStructs$est * model$Yscale / model$Xscale
      model$TreeStructs$xmin <- sapply(model$TreeStructs$xmin, function(i) {
        if (i == 0) -Inf
        else model$Xsplits[i]
      })
      model$TreeStructs$xmax <- sapply(model$TreeStructs$xmax, function(i) {
        if (i == (length(model$Xsplits) + 1)) Inf
        else model$Xsplits[i]
      })
    }

  # Fixed tree index
  } else {
    model$TreeStructs <- as.data.frame(model$TreeStructs)
    if (model$dlmType == "gp") {
      colnames(model$TreeStructs) <- c("Iter", "Tree", "fixedIdx",
                                       paste0("Lag", 1:model$pExp))
      model$TreeStructs[,-c(1:3)] <- model$TreeStructs[,-c(1:3)] *
        model$Yscale / model$Xscale
    } else {
      colnames(model$TreeStructs) <- c("Iter", "Tree", "fixedIdx", "dlnmTerm",
                                       "xmin", "xmax", "tmin", "tmax", "est")
      model$TreeStructs$est <- model$TreeStructs$est *
        model$Yscale / model$Xscale
    }
  }



  # Remove model and exposure data
  if (!save.data) {
    model$X <- NULL
    model$Tcalc <- NULL
    model$Xcalc <- NULL
    model$Z <- NULL
    model$Mo <- NULL
  }

  class(model) <- "dlmtree"
  return(model)
}
