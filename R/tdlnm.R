#' tdlnm
#'
#' @param formula object of class formula, a symbolic description of the fixed 
#' effect model to be fitted, e.g. y ~ a + b
#' @param data data frome containing variables used in the formula
#' @param exposure.data numerical matrix of complete exposure data with same
#' length as data
#' @param exposure.se numerical matrix of exposure standard errors with same
#' length as data or scalar smoothing factor representing a uniform standard
#' error over each exposure measurement
#' @param exposure.splits scalar indicating the number of splits (divided
#' evenly across quantiles of the exposure data) or list with two components:
#' 'type' set to 'values' or 'quantiles', and 'split.vals' set as a numerical
#' vector indicating the corresponding exposure values or quantiles for splits
#' @param n.trees integer for number of trees in ensemble (default = 20)
#' @param n.burn integer for length of burn-in
#' @param n.iter integer for number of iterations to run model after burn-in
#' @param n.thin integer thinning factor, i.e. keep every tenth iteration
#' @param tree.params numerical vector of alpha and beta hyperparameters 
#' controlling tree depth (see Bayesian CART, 1998)
#' @param step.prob numerical vector for probability of 1) grow/prune, 2) change
#' @param subset integer vector to analyze only a subset of data and exposures
#' @param verbose true/false print output
#' @param diagnostics true/false keep model diagnostic such as terminal nodes, 
#' acceptance details, etc.
#' @param ... NA
#'
#' @return
#' @export
#'
#' @examples
tdlnm <- function(formula, 
                  data,
                  exposure.data, 
                  exposure.se = NULL, 
                  exposure.splits = 100,
                  n.trees = 20, 
                  n.burn = 2000, 
                  n.iter = 5000, 
                  n.thin = 10,
                  tree.params = c(.95, 2), 
                  step.prob = c(.3, .4),
                  subset = c(),
                  verbose = TRUE, 
                  diagnostics = FALSE, 
                  ...)
{
  # ---- Check inputs ----
  options(stringsAsFactors = F)
  if (!is.data.frame(data)) stop("Data must be a data.frame.")
  if (!is.numeric(exposure.data)) stop("`exposure.data` must be a numeric matrix")
  if (!is.null(exposure.se)) {
    if (!is.numeric(exposure.se)) stop("`exposure.se` must be a numeric matrix")
    if (length(exposure.se) == 1)
      exposure.se <- matrix(exposure.se, nrow(exposure.data), ncol(exposure.data))
  }
  if (nrow(data) != nrow(exposure.data))
    stop("`data` and `exposure.data` must have same number of observations.")
  if (!all(sapply(list(n.trees,n.burn,n.iter,n.thin), function(i) is.numeric(i) & i>0)))
    stop("n.* must be numeric and > 0")

  if (length(tree.params) != 2)
    stop("tree.params must have length 2")
  if (any(step.prob < 0) | any(step.prob > 1))
    stop("step.prob components are incorrectly specified")
  model <- new.env(hash = T)

  # ---- Setup control and response variables ----
  model$formula <- as.formula(formula)
  tf <- terms.formula(model$formula, data = data)
  if (!attr(tf, "response")) stop("No valid response in formula.")
  if (attr(tf, "intercept")) model$intercept = TRUE
  else model$intercept = FALSE
  if (length(which(attr(tf, "term.labels") %in% colnames(data))) == 0)
    stop("No valid variables in formula.")

  # ---- Control arguments ----
  model$ctr <- list("n.trees" = n.trees, "n.burn" = n.burn,
                    "n.iter" = n.iter, "n.thin" = n.thin,
                    "verbose" = verbose, "diagnostics" = diagnostics)
  model$tree.params <- tree.params
  model$step.prob <- c(step.prob[1], step.prob[1], step.prob[2], 0)
  model$step.prob <- model$step.prob/sum(model$step.prob)
  args <- list(...)
  if (!is.null(args$par)) {
    model$ctr$par = T
    model$ctr$core = args$core
    model$ctr$rec = args$rec
  } else {
    model$ctr$par = F
    model$ctr$core = 1
    model$ctr$rec = 5000
  }

  if (model$ctr$verbose & !model$ctr$par) {
    cat("Preparing data...\n")
  }

  # ---- Create data subset ----
  if (length(subset) > 1 & is.numeric(subset)) {
    data <- data[subset,]
    exposure.data <- exposure.data[subset,]
    if (!is.null(exposure.se))
      exposure.se <- exposure.se[subset,]
    if (nrow(data) == 0)
      stop("subset did not leave any data!")
  }

  model$Mo <- dlnmMod(exposure.data, exposure.se, exposure.splits)
  data <- droplevels(data)
  model$Z <- scaleModelMatrix(model.matrix(model$formula, data))
  model$Y <- model.response(model.frame(model$formula, data = data))
  model$Y <- scale(model$Y, center = sum(range(model$Y))/2, scale = sqrt(crossprod(model$Y))[1,1])
  model$rescale$Y <- attr(model$Y, "scaled:scale")
  model$rescale$Z <- attr(model$Y, "scaled:scale") / attr(model$Z, "scaled:scale")
  rm(exposure.data, exposure.se, data)
  gc(verbose = F)


  # ---- Run model ----
  tdlnm.gaussian(model)

  # ---- Prepare output ----
  model$Z <- NULL
  rm(Z, envir = model)
  model$Y <- model$Y * model$rescale$Y + attr(model$Y, "scaled:center")

  # Change env to list
  model$Mo$Xcounts <- model$Mo$ZtX <- model$Mo$VgZtX <- NULL
  rm(Xcounts, ZtX, VgZtX, envir = model$Mo)
  Mo <- lapply(names(model$Mo), function(i) model$Mo[[i]])
  names(Mo) <- names(model$Mo)
  model$Mo <- Mo
  model.out <- lapply(names(model), function(i) model[[i]])
  names(model.out) <- names(model)

  class(model.out) <- "tdlnm"
  return(model.out)
}