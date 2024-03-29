% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transreg.R
\name{transreg}
\alias{transreg}
\title{Parametric regression models for infectious disease transmission}
\usage{
transreg(
  formula,
  sus,
  data,
  subset = NULL,
  na.action,
  dist = "weibull",
  xdist = dist,
  init = NULL,
  fixed = NULL,
  optim_method = "BFGS",
  ...
)
}
\arguments{
\item{formula}{A formula of the form "response ~ terms". The response
must be an object returned by \code{\link[survival]{Surv}}. The formula
terms can include an \code{ext} term giving an external row indicator.}

\item{sus}{A name (without quotes) or character string (with quotes)
giving the name of the variable in \code{data} that contains the
susceptible member of each pair.}

\item{data}{A data frame containing the variables named in \code{formula}.}

\item{subset}{An expression indicating which rows of \code{data} should be
included in the model fit. In stats::model.frame, subsetting is done after
evaluating the formula, which can lead to empty warnings if time == 0 or
start >= stop in data rows not included in the subset. This can be avoided
by subsetting the data in the data argument itself.}

\item{na.action}{A missing-data filter applied to \code{data} after
\code{subset}. Defaults to \code{options()$na.action}.}

\item{dist}{A string partially matching a survival time distribution from
\code{\link{transreg.distributions}} to specify the internal contact
interval distribution. Current options are \code{exponential} (the
default), \code{weibull}, and \code{loglogistic}.}

\item{xdist}{A string partially matching a survival time distribution to
specify the external contact interval distribution. By default, it
is the same as \code{dist}.}

\item{init}{A named vector of initial values for estimated coefficients.}

\item{fixed}{A named vector of fixed parameter values. These can include
terms in \code{formula}, \code{intercept}, \code{logshape} for the
internal shape parameter, and \code{xlogshape)} for the external shape
parameter.}

\item{optim_method}{The method to be used by \code{\link[stats]{optim}}.}

\item{...}{Further arguments to be passed to \code{\link[stats]{optim}}.}
}
\value{
A list with class \code{transreg} that contains the following
 objects:
 \describe{
   \item{\code{call}}{The call to \code{transreg} with complete formal
     arguments.}
   \item{\code{coefficients}}{Named vector of estimated parameters.}
   \item{\code{df}}{The number of estimated coefficients.}
   \item{\code{dist}}{String naming the internal contact interval
     distribution.}
   \item{\code{fixed}}{Named vector of fixed parameter values.}
   \item{\code{init}}{Named vector of initial parameter values.}
   \item{\code{loglik}}{The maximum log likelihood.}
   \item{\code{model_matrix}}{The data frame used to fit the model.}
   \item{\code{nlnL}}{Function for calculating the negative log
     likelihood. See details below.}
   \item{\code{optim_method}}{The method used in
     \code{\link[stats]{optim}}.}
   \item{\code{response}}{The response from \code{\link[survival]{Surv}}.}
   \item{\code{sus}}{Factor giving the susceptible member of each
     ordered pair.}
   \item{\code{var}}{The estimated variance matrix.}
   \item{\code{xdist}}{String naming the external contact interval
     distribution; \code{NULL} if model formula has no \code{ext} term.}
 }
}
\description{
Fits accelerated failure time models for infectious disease transmission
using right-censored and/or left-truncated data on contact intervals in
ordered pairs of individuals and infectious contact from external sources
with individuals.
}
\details{
\strong{Survival time distributions}
 \itemize{
   \item {The exponential distribution (\code{exponential}) has a rate
     parameter \eqn{\lambda}. Its hazard function is
     \deqn{h(t, \lambda) = \lambda.}
   }
   \item{The Weibull distribution (\code{weibull}) has a rate parameter
     \eqn{\lambda} and a shape parameter \eqn{\gamma}. Its cumulative
     hazard function is
     \deqn{H(t, \lambda, \gamma) = (\lambda t)^\gamma.}
     The exponential distribution is a Weibull distribution with
     shape \eqn{\gamma = 1}.
   }
   \item{The log-logistic distribution (\code{loglogistic}) has a rate
     parameter \eqn{lambda} and a shape parameter \eqn{\gamma}. Its
     survival function is
     \deqn{S(t, \lambda, \gamma) = \frac{1}{1 + (\lambda t)^\gamma}.}
   }
 }

 \strong{Accelerated failure time models} are log-linear models for the
 rate parameter \eqn{\lambda} in the specified survival time distribution.
 Each coefficient can be interpreted as the log rate ratio for a one-unit
 increase in the corresponding covariate, so each covariate has a
 multiplicative effect on the rate.

 \strong{Internal and external transmission models} The internal
 transmission model is for the hazard of transmission between individuals
 under observation. The external model is for hazard of transmission from
 external sources (e.g., the community or the environment) to individuals
 under observation. There are four types of covariates: internal-only,
 external-only, shared covariates with the same coefficients under both
 models, and shared covariates with (possibly) different coefficients under
 the two models. For covariates included in only one model, the covariate
 is set to zero for all data rows for the other model. Shared covariates
 with equal coefficients are reported as normal. Shared covariates that can
 have unequal coefficients are included as a shared main effect and an
 interaction with the external indicator variable. Thus, the main effect
 applies to the internal model and the main effect + interaction applies to
 the external model.
}
\references{
E Kenah (2011). Contact intervals, survival analysis of
 epidemic data, and estimation of R_0. \emph{Biostatistics} 12(3):
 548-566.
}
\author{
Eben Kenah \email{kenah.1@osu.edu}
}
