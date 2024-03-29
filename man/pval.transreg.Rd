% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transreg.R
\name{pval.transreg}
\alias{pval.transreg}
\title{p-values for estimated parameters}
\usage{
\method{pval}{transreg}(treg, parm, type = "wald")
}
\arguments{
\item{treg}{An object of class \code{transreg}.}

\item{parm}{A parameter name or vector of parameter names. If missing,
p-values are calculated for all estimated parameters.}

\item{type}{The type of p-value. Options are \code{wald} for Wald p-values
and \code{lr} for likelihood ratio p-values. The latter are more accurate
but more computationally intensive.}
}
\value{
A named vector of p-values.
}
\description{
Calculates p-values for estimated parameters from a \code{transreg} model
using a normal approximation or a likelihood ratio chi-squared statistic.
}
\author{
Eben Kenah \email{kenah.1@osu.edu}
}
