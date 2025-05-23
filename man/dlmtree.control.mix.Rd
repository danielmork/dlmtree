% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control.R
\name{dlmtree.control.mix}
\alias{dlmtree.control.mix}
\title{Control settings for dlmtree model fitting, when used for mixture models}
\usage{
dlmtree.control.mix(interactions = "noself", sparsity.prior = 1)
}
\arguments{
\item{interactions}{'noself' (default) which estimates interactions only between two 
different exposures, 'all' which also allows interactions within the same exposure, or 'none' 
which eliminates all interactions and estimates only main effects of each exposure.}

\item{sparsity.prior}{positive scalar hyperparameter for sparsity of exposures. (default: 1)}
}
\value{
list of mixture control parameters.
}
\description{
Control settings for dlmtree model fitting, when used for mixture models
}
