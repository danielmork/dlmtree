#' estDLM
#'
#' @title Calculates subgroup-specific lag effects for heterogeneous models
#' @description Method for calculating subgroup-specific lag effects for heterogeneous models: HDLM, HDLMM
#'
#' @param object an object of a model fit. Must be 'hdlm' or 'hdlmm'
#' @param new.data a data frame with new observations with the same number of modifiers
#' @param group.index a list of index (row numbers) for subgroup specification
#' @param conf.level confidence level for credible interval of effects
#' @param exposure exposure of interest for 'hdlmm' method
#' @param return.mcmc store mcmc in the output
#' @param mem.safe boolean memory parameter for rule index
#' @param verbose TRUE (default) or FALSE: print output
#'
#' @returns A list of distributed lag effects per subgroups
#' @export 
#'
estDLM <- function(object,
                   new.data,
                   group.index,
                   conf.level = 0.95,
                   exposure = NULL,
                   #cenval = 0,
                   return.mcmc = FALSE,
                   mem.safe = FALSE,
                   verbose = TRUE)
{
  if (!(class(object) %in% c("hdlm", "hdlmm", "dlmtree"))) {
    stop("`object` must be one of `hdlm` or `hdlmm`")
  }

  if (class(object) %in% c("hdlmm", "dlmtree")){
    if(!(exposure %in% object$expNames || is.null(exposure))){
      stop("exposure must match the exposure names in the model")
    } 

    if(is.null(exposure)){
      exposure <- object$expNames[1]
      warning(paste0("'exposure' has not been selected. Running with the first exposure: ", exposure))
    }

    TreeStructs <- object$TreeStructs[(object$TreeStructs$exp + 1) == which(object$expNames == exposure), ]
  } else {
    TreeStructs <- object$TreeStructs
  }
    
  if (!all(object$modNames %in% colnames(new.data))) {
    stop("`new.data` must have the same colunm names as the original model")
  }

  out             <- list()
  class(out)      <- "estDLM"
  out$conf.level  <- conf.level
  ci.lims         <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)

  # Load modifier data
  mod <- lapply(object$modNames, function(m) {
    if (!is.numeric(object$MoUnique[[m]])) {
      if (!all(unique(new.data[[m]]) %in% object$MoUnique[[m]])){
        stop("column ", m, " of `new.data` contains additional categories not in original dataset")
      }
    } else {
      if (min(new.data[[m]]) < object$MoUnique[[m]][1] || max(new.data[[m]]) > object$MoUnique[[m]][2]){
        if (verbose){
          warning("column ", m, "of `new.data` exceeds range of original dataset")
        } 
      }
    }
    new.data[[m]]
  })
  
  names(mod)  <- object$modNames
  out$mod     <- mod

  # Name group indices
  if (!is.list(group.index)) {
    group.index <- list(group.index)
  }
    
  if (is.null(names(group.index))) {
    names(group.index) <- 1:length(group.index)
  }
    
  if (any(sapply(group.index, function(x) any(x > nrow(new.data))))) {
    stop("at least one item of `group.index` exceeds range of data")
  }
    
  out$n           <- lapply(group.index, length)
  out$groupIndex  <- group.index

  tempMod         <- data.frame(Rule = unique(TreeStructs$Rule))
  rules           <- strsplit(tempMod$Rule, " & ", TRUE)
  mark            <- ceiling(length(rules) / 42)

  # Analyze trees for each group index list
  for (i in 1:length(group.index)) {
    if (names(group.index)[i] == "") {
      names(group.index[i]) <- i
    }
      
    if (verbose) {
      cat(paste0("Reanalyzing trees for group ", names(group.index)[i], "\n",
                 "% complete\n[0--------25--------50--------75--------100]\n '"))
    }
      
    modDat        <- lapply(mod, function(k) k[group.index[[i]]])
    names(modDat) <- names(mod)
    ri            <- ruleIdx(modDat, mem.safe)
    tempMod[[paste0("Weight", i)]] <- sapply(1:length(rules), function(i) {
      if (verbose && (i %% mark == 0)) {
        cat("'")
      }
      ri$returnWeight(rules[[i]])
    })
    rm(ri)

    if (verbose) {
      cat("\n")
    }
  }

  if (verbose) {
    cat("Calcuating DLMs...")
  }
    
  DLM <- merge.data.frame(TreeStructs, tempMod, by = "Rule")
  if (return.mcmc) {
    out$mcmc  <- list()
  }
    
  out$dlmMean <- list()
  out$dlmCI   <- list()
  out$dlmCum  <- list()


  # ---- dlmType: TDLM ----
  # } else if (object$dlFunction == "dlm") {
  # } else {
    for (i in 1:length(group.index)) {
      DLM$w.est <- DLM$est * DLM[[paste0("Weight", i)]]
      mcmc      <- dlmEst(as.matrix(DLM[,c("Iter", "Tree", "tmin", "tmax", "w.est")]), object$pExp, object$mcmcIter)
      if (return.mcmc) {
        out$mcmc[[names(group.index)[i]]]   <- mcmc
      }
      out$dlmMean[[names(group.index)[i]]]  <- rowMeans(mcmc)
      out$dlmCI[[names(group.index)[i]]]    <- apply(mcmc, 1, quantile, probs = ci.lims)
      out$dlmCum[[names(group.index)[i]]]   <- c(mean = mean(rowMeans(mcmc)), quantile(colSums(mcmc), probs = ci.lims))
    }
    out$dlFunction <- "dlm"

  lags          <- length(out$dlmMean[[1]])
  out$plotData  <- do.call(rbind, lapply(names(group.index), function(n) {
    data.frame(group = n, time = 1:lags, est = out$dlmMean[[n]], 
               lower = out$dlmCI[[n]][1,], upper = out$dlmCI[[n]][2,])
  }))
  
  return(out)
}


#' Calculates the weights for each modifier rule
#'
#' @description Method for calculating the weights for each modifier rule
#'
#' @param mod a list of modifier splitting rules
#' @param mem.safe boolean memory parameter
#'
#' @returns A list of weights per rule with modifiers
#'
ruleIdx <- function(mod, mem.safe = FALSE) {
  self      <- environment()
  n         <- length(mod[[1]])
  idxList   <- list()
  `%notin%` <- Negate(`%in%`)

  returnWeight = function(rules){
    if (length(rules) == 0){ return(1) }
    if (mem.safe) {
      return(length(Reduce(cppIntersection, lapply(rules, function(i) which(eval(parse(text = i))))))/self$n)
    } else {
      missingIdx  <- sapply(rules, function(i) is.null(idxList[[i]]))
      newIdx      <- list()
      idx         <- list()
      if (any(missingIdx)) {
        newIdx        <- lapply(rules[missingIdx], function(i) which(eval(parse(text = i))))
        names(newIdx) <- rules[missingIdx]
        idxList       <- append(idxList, newIdx)
      }
      return(length(Reduce(cppIntersection, lapply(rules, function(i) idxList[[i]])))/self$n)
    }
  }

  return(self)
}
