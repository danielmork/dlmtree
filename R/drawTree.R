drawTree <- function(depth, alpha, beta) {
  n.term <- 0
  lp <- log(alpha) - beta * log(1 + depth)
  if (log(runif(1)) < lp) {
    n.term <- drawTree(depth + 1, alpha, beta) +
      drawTree(depth + 1, alpha, beta)
  } else {
    n.term <- 1
  }
  return(n.term)
}
# hist(sapply(1:10000, function(i) drawTree(0, .5, .5)))
