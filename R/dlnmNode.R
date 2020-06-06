#' dlnmNode
#'
#' @param depth Depth of node in tree (starts at 0)
#' @param gen2 T/F, does node have child nodes that are both terminal
#' @param xmin Min exposure value
#' @param xmax Max exposure value
#' @param tmin Min time
#' @param tmax Max time
#' @param xsplit Split on exposure value
#' @param tsplit Split on time
#' @param counts Number of observations that are in node
#' @param logprule Log probability of split
#' @param avail.mod Numeric list of available modifiers
#' @param parent Reference to parent node
#' @param sib Reference to sibling node
#' @param c1 Reference to child node 1
#' @param c2 Reference to child node 2
#'
#' @description DLNM node environment
#'
#' @return function environment with objects and methods
#'
#'

dlnmNode <- function(depth = 0L,
                     gen2 = FALSE,
                     xmin = -Inf,
                     xmax = Inf,
                     tmin = 1L,
                     tmax = Inf,
                     xsplit = NULL,
                     tsplit = NULL,
                     counts = c(),
                     logprule = 0,
                     avail.mod = NULL,
                     parent = NULL,
                     sib = NULL,
                     c1 = NULL,
                     c2 = NULL)
{
  self <- environment()
  ZtX <- matrix()
  VgZtX <- matrix()
  update <- TRUE

  # Copy current tree structure, preserve all references
  copy = function() {
    NewTree <- dlnmNode()
    NewTree$depth <- depth
    NewTree$gen2 <- gen2
    NewTree$xmin <- xmin
    NewTree$xmax <- xmax
    NewTree$tmin <- tmin
    NewTree$tmax <- tmax
    NewTree$xsplit <- xsplit
    NewTree$tsplit <- tsplit
    NewTree$counts <- counts
    NewTree$logprule <- logprule
    NewTree$avail.mod <- avail.mod
    NewTree$ZtX <- ZtX
    NewTree$VgZtX <- VgZtX
    NewTree$update <- update
    if (!is.null(c1)) {
      NewTree$c1 <- c1$copy()
      NewTree$c2 <- c2$copy()
      NewTree$c1$sib <- NewTree$c2
      NewTree$c2$sib <- NewTree$c1
      NewTree$c1$parent <- NewTree$c2$parent <- NewTree
    }
    NewTree
  }

  # Delete all elements of tree (memory saving)
  # delete = function() {
  #   if (!is.null(c1)) {
  #     c1$delete()
  #     c2$delete()
  #     c1$parent <- c2$parent <- c1$sib <- c2$sib <- NULL
  #   }
  #   c1 <<- c2 <<- sib <<- parent <<- counts <<- avail.mod <<- self <<- NULL
  # }


  # Output references in list to all terminal nodes
  list.terminal = function() {
    t <- list()
    if (!is.null(c1)) {
      t[[1]] <- c1$list.terminal()
      t[[2]] <- c2$list.terminal()
    } else {
      t[[1]] <- self
    }
    unlist(t)
  }


  # Output references in list to all terminal nodes
  #list.terminal2 = function() {
  #  t <- list()
  #  if (!is.null(c1)) {
  #    t[[1]] <- c1
  #    t[[2]] <- c1$list.terminal2()
  #    t[[3]] <- c2$list.terminal2()
  #  }
  #  unlist(t)
  #}


  # Output ref in list to all gen2 nodes (have 2 terminal)
  list.gen2 = function() {
    t <- list()
    if (gen2) t[[1]] <- self
    if (!is.null(c1)) {
      t[[2]] <- c1$list.gen2()
      t[[3]] <- c2$list.gen2()
    }
    unlist(t)
  }


  # Output ref in list to all internal nodes (non-terminal)
  list.internal = function(min.depth = 0) {
    t <- list()
    if (!is.null(c1)) {
      if (c1$depth >= min.depth) t[[1]] <- self
      t[[2]] <- c1$list.internal(min.depth)
      t[[3]] <- c2$list.internal(min.depth)
    }
    unlist(t)
  }


  # Output ref in list to all nodes at and below current node
  list.subnodes = function() {
    t <- list()
    t[[1]] <- self
    if (!is.null(c1)) {
      t[[2]] <- c1$list.subnodes()
      t[[3]] <- c2$list.subnodes()
    }
    unlist(t)
  }


  # Output ref in list to all nodes above current node
  # list.parents = function(s = T) {
  #   if (depth == 0)
  #     return(self)
  #   else {
  #     if (s) return(unlist(parent$list.parents(F)))
  #     else return(unlist(list(self, parent$list.parents(F))))
  #   }
  # }


  # Update limits of subnodes based on observed splits
  update.subnodes = function(M) {
    counts <<- c()
    if (xmin >= xmax | tmin > tmax)
      return(FALSE)



    # Update available modifiers
    avail.mod <<- M$avail.mod(self)

    # Update subnodes
    if (!is.null(c1)) {
      c1$update <- c2$update <- TRUE
      logprule <<- M$log.pRule(self)
      if (!is.null(xsplit)) {
        c1$xmin <- xsplit; c1$xmax <- xmax
        c2$xmin <- xmin; c2$xmax <- xsplit
        c1$tmin <- c2$tmin <- tmin
        c1$tmax <- c2$tmax <- tmax
        if (!c1$update.subnodes(M)) return(FALSE)
        if (!c2$update.subnodes(M)) return(FALSE)
      } else {
        c1$tmin <- tmin; c1$tmax <- tsplit
        c2$tmin <- tsplit+1; c2$tmax <- tmax
        c1$xmin <- c2$xmin <- xmin
        c1$xmax <- c2$xmax <- xmax
        if (!c1$update.subnodes(M)) return(FALSE)
        if (!c2$update.subnodes(M)) return(FALSE)
      }
    }

    return(TRUE)
  }

  # Function to count obs in node, update if needed
  get.counts = function(model, ctr) {
    # if no update needed, return counts
    if (depth == 0)
      return(model$Mo$Xcounts[[length(model$Mo$Xcounts)]])
    if (!update)
      return(counts)

    # count exposure occurances for each node
    # use parent counts to simply update for adjacent nodes
    p.counts <- parent$get.counts(model, ctr)
    if (!is.null(c1))
      c1$update <- c2$update <- TRUE
    if (!is.null(sib$c1))
      sib$c1$update <- sib$c2$update <- TRUE

    # if parent counts update fails, then this update fails too
    if (length(p.counts) == 0)
      return(c())

    # Special case if first split is on time
    #  (then new counts are divided proportionally)
    if (xmin == -Inf && xmax == Inf) {
      if (length(model$Mo$Xsplits) == 0) {
        if (tmin == 1) {
          counts <<- model$Mo$Tsums[[tmax]][[1]]
          ZtX <<- model$Mo$ZtX[[tmax]][[1]]
          VgZtX <<- model$Mo$VgZtX[[tmax]][[1]]
        } else if (tmax == ncol(model$Mo$X)) {
          counts <<- model$Mo$Tsums[[tmin-1]][[2]]
          ZtX <<- model$Mo$ZtX[[tmin-1]][[2]]
          VgZtX <<- model$Mo$VgZtX[[tmin-1]][[2]]
        } else if (length(tmin:tmax) > 1) {
          counts <<- rowSums(model$Mo$X[,tmin:tmax])
          ZtX <<- crossprod(model$Z, counts)
          VgZtX <<- crossprod(ctr$Vg, ZtX)
        } else {
          counts <<- model$Mo$X[,tmin]
          ZtX <<- crossprod(model$Z, counts)
          VgZtX <<- crossprod(ctr$Vg, ZtX)
        }
        update <<- FALSE
      } else {
        scale <- ((tmax - tmin + 1) / ncol(model$Mo$X))
        len <- length(model$Mo$Xcounts)
        counts <<- model$Mo$Xcounts[[len]] * scale
        ZtX <<- model$Mo$ZtX[[len]] * scale
        VgZtX <<- model$Mo$VgZtX[[len]] * scale
        update <<- FALSE
      }

    # Special case if across all time
    } else if (tmin == 1 && tmax == ncol(model$Mo$X)) {
      if (xmin == -Inf) {
        id <- which(model$Mo$Xsplits == xmax)
        counts <<- model$Mo$Xcounts[[id]]
        ZtX <<- model$Mo$ZtX[[id]]
        VgZtX <<- model$Mo$VgZtX[[id]]
      } else if (xmax == Inf) {
        id <- which(model$Mo$Xsplits == xmin)
        len <- length(model$Mo$Xcounts)
        counts <<- model$Mo$Xcounts[[len]] - model$Mo$Xcounts[[id]]
        ZtX <<- model$Mo$ZtX[[len]] - model$Mo$ZtX[[id]]
        VgZtX <<- model$Mo$VgZtX[[len]] - model$Mo$VgZtX[[id]]
      } else {
        id1 <- which(model$Mo$Xsplits == xmin)
        id2 <- which(model$Mo$Xsplits == xmax)
        counts <<- model$Mo$Xcounts[[id2]] - model$Mo$Xcounts[[id1]]
        ZtX <<- model$Mo$ZtX[[id2]] - model$Mo$ZtX[[id1]]
        VgZtX <<- model$Mo$VgZtX[[id2]] - model$Mo$VgZtX[[id1]]
      }
      update <<- FALSE

    # Count manually
    } else {
      if (is.null(model$Mo$SE))
        nc <- nodeCount(model$Mo$X, model$Z, ctr$Vg, p.counts,
                        xmin, xmax, tmin, tmax)
      else
        nc <- nodeCountSE(model$Mo$X, model$Mo$SE, model$Z, ctr$Vg, p.counts,
                          xmin, xmax, tmin, tmax)

      if (!nc$Empty) {
        counts <<- nc$Count
        ZtX <<- nc$ZtX
        VgZtX <<- nc$VgZtX
        update <<- FALSE
        sib$counts <- nc$SibCount
        sib$ZtX <- nc$SibZtX
        sib$VgZtX <- nc$SibVgZtX
        sib$update <- FALSE
        return(counts)
      } else {
        return(c())
      }
    }

    # Update sibling node using parent
    if (sib$update) {
      sib$counts <- p.counts - counts
      sib$ZtX <- parent$ZtX - ZtX
      sib$VgZtX <- parent$VgZtX - VgZtX
      sib$update <- FALSE
    }

    return(counts)
  }


  # Output diagram of binary modifier tree
  # diagram = function(s = "") {
  #   if (depth == 0) cat("\nTree")
  #   if (!is.null(c1)) {
  #     cat(paste0("\n", s, "+ xs:", xsplit, ", ts:", tsplit))
  #     c1$diagram(paste0(s, "  "))
  #     c2$diagram(paste0(s, "  "))
  #   } else {
  #     cat(paste0("\n", s, "+ X:[", xmin, ", ", xmax, "] T:[", tmin, ", ", tmax, "]"))
  #   }
  # }

  # Generate p(T|alpha.new)/p(T|alpha.old)
  # pTree.alpha = function(alpha.new, alpha.old, beta, eq = "logit") {
  #   t <- list.terminal()
  #   nt <- list.internal()
  #   s <- do.call(sum, lapply(t, function(i) {
  #     logPSplit(alpha.new, beta, i$depth, TRUE, eq) -
  #       logPSplit(alpha.old, beta, i$depth, TRUE, eq)
  #   }))
  #   if (length(nt) > 0) {
  #     s <- s + do.call(sum, lapply(nt, function(i) {
  #       logPSplit(alpha.new, beta, i$depth, FALSE, eq) -
  #         logPSplit(alpha.old, beta, i$depth, FALSE, eq)
  #     }))
  #   }
  #   return(s)
  # }
  return(environment())
}
