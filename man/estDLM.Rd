% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estDLM.R
\name{estDLM}
\alias{estDLM}
\title{Calculates subgroup-specific lag effects for heterogeneous models}
\usage{
estDLM(
  object,
  new.data,
  group.index,
  conf.level = 0.95,
  exposure = NULL,
  keep.mcmc = FALSE,
  mem.safe = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{object}{object of a model fit. Must be 'hdlm' or 'hdlmm'}

\item{new.data}{data frame with new observations with the same number of modifiers}

\item{group.index}{list of index (row numbers) for subgroup specification}

\item{conf.level}{confidence level for credible interval of effects}

\item{exposure}{exposure of interest for 'hdlmm' method}

\item{keep.mcmc}{store mcmc in the output}

\item{mem.safe}{boolean memory parameter for rule index}

\item{verbose}{TRUE (default) or FALSE: print output}
}
\value{
A list with the following components:
\item{conf.level}{Specified confidence level}
\item{mod}{a list of modifers with a vector of values from the model}
\item{n}{Number of observation per specified subgroup}
\item{groupIndex}{list of index (row numbers) for specified subgroup}
\item{dlmMean}{distributed lag effects per subgroups}
\item{dlmCI}{credible intervals for distributed lag effects per subgroups}
\item{dlmCum}{cumulative effects per subgroups}
\item{dlFunction}{type of DLM class}
\item{plotData}{data frame built for easier visualization of distributed lag effects for each subgroup (facet)}
}
\description{
Method for calculating subgroup-specific lag effects for heterogeneous models: HDLM, HDLMM
}
\details{
estDLM
}
