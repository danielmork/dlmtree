#' Posterior inclusion probabilities (PIPs) for modifiers in HDLM
#'
#' @param object An object of class dlmtree.
#' @param type Type=1 indicates single modifier PIPs. Type=2 indicates joint modifier PIPs for two modifiers.
#'
#' @return A vector (type=1) or data.frame (type=2) of PIPs.
#' @export


pip <- function(object, type=1){
  
  
  if(type==1){ # main effect PIPs
    
    return(colMeans(object$modCount>0))
    
  }else if(type==2){ # interaction PIPs
    
    
    sp <- cbind.data.frame(Rule = object$termRules, object$TreeStructs[,2:4])
    sp <- sp[!duplicated(sp),]
    splitRules <- lapply(strsplit(sp$Rule, "&", T), function(i) {
      sort(as.numeric(sapply(strsplit(i, ">=|<|\\[\\]|\\]\\[", perl = T), function(j) j[1])))
    })
    splitCount <- lapply(1:object$mcmcIter, function(i) list())
    treeMods <- lapply(1:object$mcmcIter, function(i) rep(0, object$nTrees))
    for (i in 1:length(splitRules)) {
      it <- sp$Iter[i]
      tr <- sp$Tree[i]
      n <- length(splitRules[[i]])
      treeMods[[it]][tr+1] <- max(treeMods[[it]][tr+1], length(unique(splitRules[[i]])))
      if (n == 1) {
      } else if (n > 1) {
        for (s in 1:(n-1)) {
          for (e in (s+1):n) {
            splitCount[[it]][[paste0(object$modNames[splitRules[[i]][s]+1], ".", 
                                     object$modNames[splitRules[[i]][e]+1])]] <- T
          }
        }
      }
    }
    mean(do.call(rbind, treeMods)==4)
    sc <- do.call(bind_rows, lapply(splitCount, cbind.data.frame))
    sc.mean <- sort(colMeans(!is.na(sc)))
    sc.mat <- data.frame()
    for (i in 1:length(sc.mean)) {
      names <- sort(strsplit(names(sc.mean)[i], ".", T)[[1]])
      sc.mat <- rbind.data.frame(sc.mat,
                                 data.frame("var1" = names[1], "var2" = names[2], "pip" = sc.mean[i]),
                                 data.frame("var1" = names[2], "var2" = names[1], "pip" = sc.mean[i]))
    }
    
    # remove duplciates and row names
    rownames(sc.mat) <- NULL
    sc.mat <- unique(sc.mat)
    
    # this adds a zero for any combinations that did not appear
    sc.mat <- merge(expand.grid(var1=object$modNames, var2=object$modNames),sc.mat, by=c("var1","var2"), all=TRUE)
    if(any(is.na(sc.mat$pip))){
      sc.mat$pip[which(is.na(sc.mat$pip))] <- 0
    }
    
    # sort from largest to smallest PIP
    sc.mat <- sc.mat[order(-sc.mat$pip),]
    
    return(sc.mat)
  }
  
}




