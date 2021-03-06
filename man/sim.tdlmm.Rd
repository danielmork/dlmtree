% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.tdlmm.R
\name{sim.tdlmm}
\alias{sim.tdlmm}
\title{sim.tdlmm}
\usage{
sim.tdlmm(
  sim = 1,
  error = 10,
  mean.p = 0.5,
  n.exp = 25,
  prop.active = 0.05,
  n = 5000,
  expList = NULL
)
}
\arguments{
\item{sim}{integer (1-4) specifying simulation scenario}

\item{error}{positive scalar specifying error variance}

\item{mean.p}{scalar between zero and one specifying mean probability
for simulation scenario one}

\item{n.exp}{number of exposures for simulation scenarios three and four}

\item{prop.active}{proportion of active exposures for simulation scenario
three}

\item{n}{sample size for simulation}

\item{expList}{named list of exposure data}
}
\value{

}
\description{
Simulation scenarios to accompany TDLM/TDLMM
}
