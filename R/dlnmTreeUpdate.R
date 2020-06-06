#' dlnmTreeUpdate
#'
#' @param tree Reference to tree environment
#' @param step Step (1-4 relates to grow, prune, change, swap)
#' @param model Model environment
#'
#' @description Propose an update for dlnmTree
#'
#' @return List detailing update
#'
dlnmTreeUpdate <- function(tree, step, model)
{
  term <- tree$list.terminal()
  ret <- list("Accept" = 0, "Step" = "", "MHR" = -Inf, "MHR2" = -Inf,
              "Covar" = "", "Val" = 0, "Depth" = 0, "nTerm" = length(term))
  # ---- Grow two new terminal nodes ----
  if (step == 1 | length(term) == 1) {
    ret$Step = "Grow"
    # Select terminal node
    node <- term[[sample.int(length(term), 1)]]
    ret$Depth <- node$depth

    # Propose split
    new.split <- model$Mo$propose.split(node)
    if (!new.split$Success) return(ret)

    # Add child nodes
    node$logprule <- new.split$logprule
    if (new.split$mod == "X") {
      node$xsplit <- new.split$split
      node$c1 <- dlnmNode(node$depth + 1L, FALSE, new.split$split, node$xmax,
                          node$tmin, node$tmax, NULL, NULL, c(), 0, NULL, node)
      node$c2 <- dlnmNode(node$depth + 1L, FALSE, node$xmin, new.split$split,
                          node$tmin, node$tmax, NULL, NULL, c(), 0, NULL, node)
    } else {
      node$tsplit <- new.split$split
      node$c1 <- dlnmNode(node$depth + 1L, FALSE, node$xmin, node$xmax,
                          node$tmin, new.split$split, NULL, NULL, c(), 0, NULL, node)
      node$c2 <- dlnmNode(node$depth + 1L, FALSE, node$xmin, node$xmax,
                          new.split$split + 1, node$tmax, NULL, NULL, c(), 0, NULL, node)
    }
    node$c2$sib <- node$c1
    node$c1$sib <- node$c2
    node$gen2 <- TRUE
    if (node$depth > 0)
      node$parent$gen2 <- FALSE
    if (!node$c1$update.subnodes(model$Mo) | !node$c2$update.subnodes(model$Mo))
      return(ret)

    # Calculate p(grow)
    ret$MHR <- log(length(term)) - log(length(tree$list.gen2())) +
      2 * logPSplit(model$tree.params[1], model$tree.params[2],
                    node$depth + 1, TRUE) +
      logPSplit(model$tree.params[1], model$tree.params[2], node$depth, FALSE) -
      logPSplit(model$tree.params[1], model$tree.params[2], node$depth, TRUE)

    ret$Accept = 1
    ret$Covar = new.split$mod
    ret$Val = new.split$split
    return(ret)


  # ---- Prune two terminal nodes ----
  } else if (step == 2) {
    ret$Step = "Prune"
    # Select gen2 node to prune
    gen2 <- tree$list.gen2()
    node <- gen2[[sample.int(length(gen2), 1)]]
    ret$Depth <- node$depth

    # Calculate MHR
    ret$MHR <- log(length(gen2)) - log(length(term) - 1) -
      2 * logPSplit(model$tree.params[1], model$tree.params[2],
                    node$depth + 1, TRUE) -
      logPSplit(model$tree.params[1], model$tree.params[2], node$depth, FALSE) +
      logPSplit(model$tree.params[1], model$tree.params[2], node$depth, TRUE)

    # Get split details
    if (!is.null(node$xsplit)) {
      ret$Covar <- "X"
      ret$Val <- node$xsplit
    } else {
      ret$Covar <- "T"
      ret$Val <- node$tsplit
    }

    # Prune node
    node$c1$counts <- node$c2$counts <- c()
    node$c1 <- node$c2 <- NULL
    node$logprule <- 0
    node$gen2 <- F
    node$xsplit <- node$tsplit <- NULL
    if (is.null(node$sib$c1) && node$depth > 0)
      node$parent$gen2 <- TRUE

    ret$Accept = 1
    return(ret)


  # ---- Change rule at internal node ----
  } else if (step == 3) {
    ret$Step = "Change"
    # Select node to change
    nt <- tree$list.internal()
    node <- nt[[sample.int(length(nt), 1)]]
    ret$Depth <- node$depth
    sub.nodes <- node$list.subnodes()

    # Calculate old p(rule)
    ret$MHR <- node$logprule - do.call(sum, lapply(sub.nodes, function(i) i$logprule))

    # Generate new split
    new.split <- model$Mo$propose.split(node)
    if (!new.split$Success) return(ret)

    # Update node
    node$logprule <- new.split$logprule
    node$c1$update <- node$c2$update <- TRUE
    ret$Covar <- new.split$mod
    ret$Val <- new.split$split
    if (new.split$mod == "X") {
      ret$Covar <- "X"
      node$xsplit <- ret$Val <- new.split$split
      node$tsplit <- NULL
      node$c1$xmin <- node$xsplit
      node$c1$xmax <- node$xmax
      node$c2$xmin <- node$xmin
      node$c2$xmax <- node$xsplit
      node$c1$tmin <- node$c2$tmin <- node$tmin
      node$c1$tmax <- node$c2$tmax <- node$tmax
      if (!node$c1$update.subnodes(model$Mo) |
          !node$c2$update.subnodes(model$Mo))
        return(ret)
    } else {
      ret$Covar <- "T"
      node$tsplit <- ret$Val <- new.split$split
      node$xsplit <- NULL
      node$c1$tmin <- node$tmin
      node$c1$tmax <- node$tsplit
      node$c2$tmin <- node$tsplit+1
      node$c2$tmax <- node$tmax
      node$c1$xmin <- node$c2$xmin <- node$xmin
      node$c1$xmax <- node$c2$xmax <- node$xmax
      if (!node$c1$update.subnodes(model$Mo) |
          !node$c2$update.subnodes(model$Mo))
        return(ret)
    }
    # Calculate new p(rule)
    ret$MHR <- ret$MHR - node$logprule + do.call(sum, lapply(sub.nodes, function(i) i$logprule))
    ret$Accept = 1
    return(ret)



  } else {
    ret$Step = "Shift"

    # Find all internal nodes
    nt <- tree$list.internal()
    if (length(nt) == 0) return(ret)
    node <- nt[[sample.int(length(nt), 1)]]
    node$c1$update <- node$c2$update <- TRUE
    sub.nodes <- node$list.subnodes()
    ret$MHR <- -do.call(sum, lapply(sub.nodes, function(i) i$logprule))

    # Perform T shift
    if (!is.null(node$tsplit)) {
      ret$Covar <- "T"
      old.split <- node$tsplit
      if (old.split == 1) {
        new.split <- min(2, ncol(model$Mo$X) - 1)
        if (new.split < ncol(model$Mo$X) - 1)
          ret$MHR <- ret$MHR - log(2)
      } else if (old.split == (ncol(model$Mo$X) - 1)) {
        new.split <- max(1, old.split - 1)
        if (new.split > 1)
          ret$MHR <- ret$MHR - log(2)
      } else {
        new.split <- old.split + (-1)^rbinom(1, 1, .5)
        if (new.split == 1 || new.split == (ncol(model$Mo$X) - 1))
          ret$MHR <- ret$MHR + log(2)
      }
      node$tsplit <- new.split
      ret$Val <- new.split

    # Perform X shift
    } else {
      ret$Covar <- "X"
      old.split <- node$xsplit
      old.split.idx <- which(model$Mo$Xsplits == old.split)
      max.xsplit <- max(model$Mo$Xsplits)
      if (old.split == model$Mo$Xsplits[1]) {
        new.split <- min(model$Mo$Xsplits[2], max.xsplit)
        if (new.split < max.xsplit)
          ret$MHR <- ret$MHR - log(2)
      } else if (old.split == max.xsplit) {
        new.split.idx <- max(1, old.split.idx - 1)
        new.split <- model$Mo$Xsplits[new.split.idx]
        if (new.split.idx > 1)
          ret$MHR <- ret$MHR - log(2)
      } else {
        new.split.idx <- old.split.idx + (-1)^rbinom(1, 1, .5)
        new.split <- model$Mo$Xsplits[new.split.idx]
        if (new.split.idx == 1 || new.split == max.xsplit)
          ret$MHR <- ret$MHR + log(2)
        node$xsplit <- new.split
        ret$Val <- new.split
      }
    }

    if (!node$update.subnodes(model$Mo))
      return(ret)
    # Calculate new p(rule)
    ret$MHR <- ret$MHR + do.call(sum, lapply(sub.nodes, function(i) i$logprule))
    ret$Accept = 1
    return(ret)
  }

  # ---- Swap two internal nodes -----
  # } else {
  #   ret$Step = "Swap"
  #   if (length(term) < 3) return(ret)
  #
  # # Find all internal nodes
  # nt <- tree$list.internal(2)
  # if (length(nt) == 0) return(ret)
  # node <- nt[[sample.int(length(nt), 1)]]
  # sub.nodes <- node$parent$list.subnodes()
  # ret$MHR <- -do.call(sum, lapply(sub.nodes, function(i) i$logprule))
  #
  #   # Get current states
  #   p.xsplit <- node$parent$xsplit; p.tsplit <- node$parent$tsplit
  #   c.xsplit <- node$xsplit; c.tsplit <- node$tsplit
  #
  #   # Swap states
  #   node$parent$xsplit <- c.xsplit; node$parent$tsplit <- c.tsplit
  #   node$xsplit <- p.xsplit; node$tsplit <- p.tsplit
  #
  #   # Swap children of sibling if same split
  #   if (!is.null(node$sib$c1)) {
  #     if (!is.null(node$sib$xsplit) && !is.null(c.xsplit)) {
  #       if (node$sib$xsplit == c.xsplit) {
  #         node$sib$xsplit <- p.xsplit; node$sib$tsplit <- p.tsplit
  #       }
  #     } else if (!is.null(node$sib$tsplit) && !is.null(c.tsplit)) {
  #       if (node$sib$tsplit == c.tsplit) {
  #         node$sib$xsplit <- p.xsplit; node$sib$tsplit <- p.tsplit
  #       }
  #     }
  #   }
  #   if (!node$parent$update.subnodes(model$Mo)) return(ret)
  #
  #   ret$MHR <- ret$MHR + do.call(sum, lapply(sub.nodes, function(i) i$logprule))
  #   ret$Accept = 1
  #   return(ret)
  # }
}
