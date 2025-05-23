% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R, R/print.dlmtree.R,
%   R/print.summary.hdlm.R, R/print.summary.hdlmm.R, R/print.summary.monotone.R,
%   R/print.summary.tdlm.R, R/print.summary.tdlmm.R, R/print.summary.tdlnm.R
\name{print}
\alias{print}
\alias{print.tdlnm}
\alias{print.tdlm}
\alias{print.tdlmm}
\alias{print.hdlm}
\alias{print.hdlmm}
\alias{print.monotone}
\alias{print.summary.hdlm}
\alias{print.summary.hdlmm}
\alias{print.summary.monotone}
\alias{print.summary.tdlm}
\alias{print.summary.tdlmm}
\alias{print.summary.tdlnm}
\title{print}
\usage{
print(x, ...)

\method{print}{tdlnm}(x, ...)

\method{print}{tdlm}(x, ...)

\method{print}{tdlmm}(x, ...)

\method{print}{hdlm}(x, ...)

\method{print}{hdlmm}(x, ...)

\method{print}{monotone}(x, ...)

\method{print}{summary.hdlm}(x, digits = 3, ...)

\method{print}{summary.hdlmm}(x, digits = 3, ...)

\method{print}{summary.monotone}(x, digits = 3, ...)

\method{print}{summary.tdlm}(x, digits = 3, ...)

\method{print}{summary.tdlmm}(x, digits = 3, ...)

\method{print}{summary.tdlnm}(x, digits = 3, ...)
}
\arguments{
\item{x}{An object of class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone', 
representing a fitted model using dlmtree(); or a summary object produced by applying 
summary() to one of these model objects.}

\item{...}{additional parameters}

\item{digits}{number of decimal places to round the numeric values to}
}
\value{
For a fitted model object, prints an assorted model output including model formula call and available methods.
For a summary object, prints a summary output of a model fit in the R console.
}
\description{
print generic function for S3method
}
