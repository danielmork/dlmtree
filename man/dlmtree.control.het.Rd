% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control.R
\name{dlmtree.control.het}
\alias{dlmtree.control.het}
\title{Control settings for dlmtree model fitting, when used for heterogeneous models}
\usage{
dlmtree.control.het(
  modifiers = "all",
  modifier.splits = 20,
  modtree.params = c(0.95, 2),
  modtree.step.prob = c(0.25, 0.25, 0.25),
  dlmtree.type = "shared",
  selection.prior = 0.5
)
}
\arguments{
\item{modifiers}{string vector containing desired modifiers to be included in a modifier tree.
The strings in the vector must match the names of the columns of the data. 
By default, a modifier tree considers all covariates in the formula as modifiers unless stated otherwise.}

\item{modifier.splits}{integer value to determine the possible number of splitting points that will be used for a modifier tree.}

\item{modtree.params}{numerical vector of alpha and beta hyperparameters
controlling modifier tree depth. (default: alpha = 0.95, beta = 2)}

\item{modtree.step.prob}{numerical vector for probability of each step for modifier tree updates: 1) grow, 2) prune,
3) change. (default: c(0.25, 0.25, 0.25))}

\item{dlmtree.type}{specification of dlmtree type for HDLM: shared (default) or nested.}

\item{selection.prior}{scalar hyperparameter for sparsity of modifiers. Must be between 0.5 and 1. 
Smaller value corresponds to increased sparsity of modifiers.}
}
\value{
list of control parameters for heterogeneous models.
}
\description{
Control settings for dlmtree model fitting, when used for heterogeneous models
}
