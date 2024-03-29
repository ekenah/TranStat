#' Parametric regression models for infectious disease transmission
#'
#' Fits accelerated failure time models for infectious disease transmission
#' using right-censored and/or left-truncated data on contact intervals in
#' ordered pairs of individuals and infectious contact from external sources
#' with individuals.
#'
#' @param formula A formula of the form "response ~ terms". The response
#'  must be an object returned by \code{\link[survival]{Surv}}. The formula
#'  terms can include an \code{ext} term giving an external row indicator.
#' @param sus A name (without quotes) or character string (with quotes)
#'  giving the name of the variable in \code{data} that contains the
#'  susceptible member of each pair.
#' @param data A data frame containing the variables named in \code{formula}.
#' @param subset An expression indicating which rows of \code{data} should be
#'  included in the model fit. In stats::model.frame, subsetting is done after
#'  evaluating the formula, which can lead to empty warnings if time == 0 or
#'  start >= stop in data rows not included in the subset. This can be avoided
#'  by subsetting the data in the data argument itself.
#' @param na.action A missing-data filter applied to \code{data} after
#'  \code{subset}. Defaults to \code{options()$na.action}.
#' @param dist A string partially matching a survival time distribution from
#'  \code{\link{transreg.distributions}} to specify the internal contact
#'  interval distribution. Current options are \code{exponential} (the
#'  default), \code{weibull}, and \code{loglogistic}.
#' @param xdist A string partially matching a survival time distribution to
#'  specify the external contact interval distribution. By default, it
#'  is the same as \code{dist}.
#' @param init A named vector of initial values for estimated coefficients.
#' @param fixed A named vector of fixed parameter values. These can include
#'  terms in \code{formula}, \code{intercept}, \code{logshape} for the
#'  internal shape parameter, and \code{xlogshape)} for the external shape
#'  parameter.
#' @param optim_method The method to be used by \code{\link[stats]{optim}}.
#' @param ... Further arguments to be passed to \code{\link[stats]{optim}}.
#'
#' @details
#'  \strong{Survival time distributions}
#'  \itemize{
#'    \item {The exponential distribution (\code{exponential}) has a rate
#'      parameter \eqn{\lambda}. Its hazard function is
#'      \deqn{h(t, \lambda) = \lambda.}
#'    }
#'    \item{The Weibull distribution (\code{weibull}) has a rate parameter
#'      \eqn{\lambda} and a shape parameter \eqn{\gamma}. Its cumulative
#'      hazard function is
#'      \deqn{H(t, \lambda, \gamma) = (\lambda t)^\gamma.}
#'      The exponential distribution is a Weibull distribution with
#'      shape \eqn{\gamma = 1}.
#'    }
#'    \item{The log-logistic distribution (\code{loglogistic}) has a rate
#'      parameter \eqn{lambda} and a shape parameter \eqn{\gamma}. Its
#'      survival function is
#'      \deqn{S(t, \lambda, \gamma) = \frac{1}{1 + (\lambda t)^\gamma}.}
#'    }
#'  }
#'
#'  \strong{Accelerated failure time models} are log-linear models for the
#'  rate parameter \eqn{\lambda} in the specified survival time distribution.
#'  Each coefficient can be interpreted as the log rate ratio for a one-unit
#'  increase in the corresponding covariate, so each covariate has a
#'  multiplicative effect on the rate.
#'
#'  \strong{Internal and external transmission models} The internal
#'  transmission model is for the hazard of transmission between individuals
#'  under observation. The external model is for hazard of transmission from
#'  external sources (e.g., the community or the environment) to individuals
#'  under observation. There are four types of covariates: internal-only,
#'  external-only, shared covariates with the same coefficients under both
#'  models, and shared covariates with (possibly) different coefficients under
#'  the two models. For covariates included in only one model, the covariate
#'  is set to zero for all data rows for the other model. Shared covariates
#'  with equal coefficients are reported as normal. Shared covariates that can
#'  have unequal coefficients are included as a shared main effect and an
#'  interaction with the external indicator variable. Thus, the main effect
#'  applies to the internal model and the main effect + interaction applies to
#'  the external model.
#'
#' @return A list with class \code{transreg} that contains the following
#'  objects:
#'  \describe{
#'    \item{\code{call}}{The call to \code{transreg} with complete formal
#'      arguments.}
#'    \item{\code{coefficients}}{Named vector of estimated parameters.}
#'    \item{\code{df}}{The number of estimated coefficients.}
#'    \item{\code{dist}}{String naming the internal contact interval
#'      distribution.}
#'    \item{\code{fixed}}{Named vector of fixed parameter values.}
#'    \item{\code{init}}{Named vector of initial parameter values.}
#'    \item{\code{loglik}}{The maximum log likelihood.}
#'    \item{\code{model_matrix}}{The data frame used to fit the model.}
#'    \item{\code{nlnL}}{Function for calculating the negative log
#'      likelihood. See details below.}
#'    \item{\code{optim_method}}{The method used in
#'      \code{\link[stats]{optim}}.}
#'    \item{\code{response}}{The response from \code{\link[survival]{Surv}}.}
#'    \item{\code{sus}}{Factor giving the susceptible member of each
#'      ordered pair.}
#'    \item{\code{var}}{The estimated variance matrix.}
#'    \item{\code{xdist}}{String naming the external contact interval
#'      distribution; \code{NULL} if model formula has no \code{ext} term.}
#'  }
#'
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @references E Kenah (2011). Contact intervals, survival analysis of
#'  epidemic data, and estimation of R_0. \emph{Biostatistics} 12(3):
#'  548-566.
#' @export
transreg <- function(formula, sus, data, subset=NULL, na.action,
                     dist="weibull", xdist = dist, init=NULL, fixed=NULL,
                     optim_method="BFGS", ...) {
  # fit accelerated failure time model using pairwise data

  # match arguments and ensure that formula argument is provided
  mcall <- match.call(expand.dots = FALSE)
  indx <- match(c("formula", "data", "subset", "na.action"),
                names(mcall), nomatch = 0)
  if (indx[1] == 0) stop("A formula argument is required.")

  # pass model frame arguments to stats::model.frame
  model <- mcall[c(1, indx)]
  model[[1]] <- quote(stats::model.frame)
  special <- c("ext")
  model$formula <- terms(formula, special)
  mframe <- eval.parent(model)

  # get model responses and model matrix
  mterms <- attr(mframe, "terms")
  x <- model.matrix(mterms, mframe)
  colnames(x)[colnames(x) == "(Intercept)"] <- "intercept"
  y <- model.response(mframe, "numeric")
  ymat <- data.frame(as.matrix(y))
  attr(ymat, "type") <- attr(y, "type")

  # indicator variable for external rows
  if (is.null(attr(mterms, "specials")$ext)) {
    ext <- rep(0, nrow(x))
    xdist <- NULL
  } else {
    extname <- survival::untangle.specials(mterms, "ext")$vars
    ext <- x[, extname]
    if (all(ext == 0)) {
      stop("Formula includes ext() term, but data has no external rows.")
    } else {
      x[, "intercept"] <- ifelse(ext, 0, 1)
      colnames(x)[colnames(x) == extname] <- "xintercept"
    }
  }
  ymat$ext <- ext

  # identify susceptibles in pairs with possible transmission
  if (missing(sus)) stop("Susceptible identifier not specified.")
  if (class(substitute(sus)) == "character") sus <- as.name(sus)
  sus <- eval(substitute(sus), data)
  if (missing(subset)) {
    ymat$sus <- sus
  } else {
    ymat$sus <- sus[eval(substitute(subset), data)]
  }

  # find distribution from partial match
  dist_table <- c("exponential", "loglogistic", "weibull")
  dist <- dist_table[pmatch(dist, dist_table)]
  if (is.na(dist)) stop("Failure time distribution not recognized.")

  # initial coefficient vector with log shape parameters as needed
  beta <- rep(0, ncol(x))
  names(beta) <- colnames(x)
  if (dist == "exponential") {
    pvec <- beta
  } else {
    pvec <- c("logshape" = 0, beta)
  }
  if (!is.null(xdist) && xdist != "exponential") {
    pvec <- c(pvec, "xlogshape" = 0)
  }

  # process user-specified initial values
  if (!is.null(init)) {
    # check validity of initial value vector
    if (is.null(names(init))) stop("Initial values must be a named vector.")
    if (any(duplicated(names(init)))) {
      stop("Duplicate names in initial value vector.")
    }
    init_indices <- match(names(init), names(pvec), nomatch = 0)
    if (any(init_indices) == 0) {
      stop("Each initial value name must match a coefficient name.")
    }
    pvec[init_indices] <- init
  }

  # process user-specified fixed values
  if (!is.null(fixed)) {
    # check validity of fixed value vector
    if (is.null(names(fixed))) stop("Fixed values must be a named vector.")
    if (any(duplicated(names(fixed)))) {
      stop("Duplicate names in fixed value vector.")
    }

    # separate fixed coefficients from coefficients to be estimated
    fixed_indices <- match(names(fixed), names(pvec), nomatch = 0)
    if (any(fixed_indices == 0)) {
      stop("Each fixed value name must match a coefficient name.")
    }
    pvec <- pvec[-fixed_indices]
  }

  # full model negative log likelihood
  # fvec is given as a formal argument to allow passage through stats::optim
  nlnL <- function(pvec, fvec=fixed) {
    beta <- c(pvec, fvec)

    # get internal log shape parameter and remove it from beta
    beta_lsindex <- match("logshape", names(beta), nomatch = 0)
    if (beta_lsindex > 0) {
      lnshape <- beta["logshape"]
      beta <- beta[-beta_lsindex]
    } else {
      lnshape <- 0
    }

    # get external log shape parameter and remove it from beta
    beta_xlsindex <- match("xlogshape", names(beta), nomatch = 0)
    if (beta_xlsindex > 0) {
      lnxshape <- beta["xlogshape"]
      beta <- beta[-beta_xlsindex]
    } else {
      lnxshape <- 0
    }

    # add rate and shape parameters to ymat
    ymat$lnshape <- ifelse(ext, lnxshape, lnshape)
    ymat$lnrate <- x[, names(beta), drop = FALSE] %*% beta

    return(transreg.nlnL(ymat, dist, xdist))
  }

  # full model fit with point estimates and variance matrix
  fit <- stats::optim(pvec, nlnL, fvec = fixed, method = optim_method,
                      hessian = TRUE, ...)
  if (length(fit$hessian) > 0) {
    coef <- fit$par
    var <- solve(fit$hessian)
  } else {
    coef <- NULL
    var <- NULL
  }

  output <- structure(list(
                           call = mcall,
                           coefficients = coef,
                           df = length(coef),
                           dist = dist,
                           fixed = fixed,
                           init = init,
                           loglik = -fit$value,
                           model_matrix = x,
                           nlnL = nlnL,
                           optim_method = optim_method,
                           response = y,
                           sus = sus,
                           var = var,
                           xdist = xdist),
                      class = "transreg")
  return(output)
}

# Internal methods used by transreg ---------------------------------------
# listed in alphabetical order

transreg.distributions <- list(
  exponential = list(lnhaz = function(t, lnrate, lnshape) lnrate,
                     cumhaz = function(t, lnrate, lnshape) exp(lnrate) * t),
  loglogistic = list(lnhaz = function(t, lnrate, lnshape) {
                       rate <- exp(lnrate)
                       shape <- exp(lnshape)
                       (lnshape + shape * lnrate + (shape - 1) * log(t)
                        - log(1 + (rate * t)^shape))
                     },
                     cumhaz = function(t, lnrate, lnshape) {
                       log(1 + (exp(lnrate) * t)^exp(lnshape))
                     }),
  weibull = list(lnhaz = function(t, lnrate, lnshape) {
                   shape <- exp(lnshape)
                   lnshape + shape * lnrate + (shape - 1) * log(t)
                 },
                 cumhaz = function(t, lnrate, lnshape) {
                   (exp(lnrate) * t)^exp(lnshape)
                 })
)

transreg.nlnL <- function(ymat, dist, xdist) {
  # internal hazard and cumulative hazard functions
  lnhaz <- transreg.distributions[[dist]]$lnhaz
  cumhaz <- transreg.distributions[[dist]]$cumhaz

  # external hazard and cumulative hazard functions
  if (!is.null(xdist)) {
    xlnhaz <- transreg.distributions[[xdist]]$lnhaz
    xcumhaz <- transreg.distributions[[xdist]]$cumhaz
  }

  # get follow-up times
  hmat <- subset(ymat, status == 1)
  if (attr(ymat, "type") == "right") {
    htimes <- hmat$time
    stimes <- ymat$time
  } else if (attr(ymat, "type") == "counting") {
    htimes <- hmat$stop
    stimes <- ymat$stop
  }

  # calculate log hazards
  lnhazards <- lnhaz(t = htimes, lnrate = hmat$lnrate, lnshape = hmat$lnshape)
  if (!is.null(xdist)) {
    lnhazards <- ifelse(hmat$ext,
      xlnhaz(t = htimes, lnrate = hmat$lnrate, lnshape = hmat$lnshape),
      lnhazards
    )
  }

  # calculate cumulative hazards
  cumhazards <- cumhaz(t = stimes, lnrate = ymat$lnrate,
                       lnshape = ymat$lnshape)
  if (!is.null(xdist)) {
    cumhazards <- ifelse(ymat$ext,
                         xcumhaz(t = stimes, lnrate = ymat$lnrate,
                                 lnshape = ymat$lnshape),
                         cumhazards)
  }
  if (attr(ymat, "type") == "counting") {
    initcumhaz <- cumhaz(t = ymat$start, lnrate = ymat$lnrate,
                         lnshape = ymat$lnshape)
    if (!is.null(xdist)) {
      initcumhaz <- ifelse(ymat$ext,
                           xcumhaz(t = ymat$start, lnrate = ymat$lnrate,
                                   lnshape = ymat$lnshape),
                           initcumhaz)
    }
    cumhazards <- cumhazards - initcumhaz
  }

  lnhsums <- tapply(lnhazards, hmat$sus, matrixStats::logSumExp)
  lnh <- sum(lnhsums)
  lnS <- -sum(cumhazards)
  return(-lnh - lnS)
}

# Methods for transreg objects ============================================
# listed in alphabetical order

#' Confidence intervals for estimated parameters
#'
#' Calculates confidence intervals for estimated parameters from a
#' \code{transreg} model by inverting the Wald or likelihood ratio tests.
#'
#' @param treg An object of class \code{transeg}.
#' @param parm A parameter name or vector of parameter names. If missing,
#'  confidence intervals are calculated for all estimated parameters.
#' @param level The confidence level (1 - \eqn{alpha}).
#' @param type The type of confidence interval. Current options are \code{wald}
#'  for Wald confidence limits and \code{lr} for likelihood ratio confidence
#'  limits. The latter are more accurate but more computationally intensive.
#' @param lrint The highest and lowest multiples of the standard error used by
#'  \code{\link[stats]{uniroot}} to find likelihood ratio confidence limits.
#' @param extendInt_lower The extendInt argument in \code{\link[stats]{uniroot}}
#'  when calculating the lower likelihood ratio confidence limit.
#' @param extendInt_upper The extendInt argument in \code{\link[stats]{uniroot}}
#'  when calculating the upper likelihood ratio confidence limit.
#' @param ... Additional arguments sent to \code{\link[stats]{uniroot}}.
#'
#' @return A data frame with one row for each parameter in \code{parm} that
#'  contains the lower and upper confidence limits in columns labeled
#'  \eqn{\frac{\alpha}{2}} and \eqn{1 - \frac{\alpha}{2}} expressed as
#'  percentiles.
#'
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
confint.transreg <- function(treg, parm, level=0.95, type="wald",
                             lrint=c(.2, 5), extendInt_lower = "downX",
                             extendInt_upper = "upX", ...) {
  # validate parameters
  if (missing(parm)) {
    parm <- names(treg$coefficients)
  } else if (anyNA(match(parm, names(treg$coefficients)))) {
    stop("Each parameter must match a coefficient name.")
  }

  # get percentages for lower and upper limits
  alpha <- 1 - level
  prb <- c(alpha / 2, 1 - alpha / 2)
  pct <- paste(signif(100 * prb, digits = 3), "%", sep = "")

  # determine type
  type_table <- c("wald", "lr")
  type <- type_table[pmatch(type, type_table)]
  if (is.na(type)) stop("Type not recognized.")

  # find Z scores and standard error
  z <- qnorm(1 - alpha / 2)
  se <- sqrt(diag(treg$var)[parm])

  # confidence limits
  if (type == "wald") {
    # wald confidence limits
    lower <- treg$coefficients[parm] - z * se
    upper <- treg$coefficients[parm] + z * se
  } else {
    # likelihood ratio confidence limits
    d <- qchisq(level, df = 1)
    limits <- function(parm) {
      pvec <- treg$coefficients[-match(parm, names(treg$coefficients))]
      parm_d <- function(val) {
        fixed <- c(treg$fixed, val)
        names(fixed) <- c(names(treg$fixed), parm)
        parm_fit <- stats::optim(pvec, treg$nlnL, fvec = fixed,
                                 method = treg$optim_method)
        return(2 * (treg$loglik + parm_fit$val) - d)
      }
      lower <- list(root = -Inf)
      try(lower <- uniroot(parm_d,
                           treg$coefficients[parm] - lrint * z * se[parm],
                           extendInt = extendInt_lower, ...))
      upper <- list(root = Inf)
      try(upper <- uniroot(parm_d,
                           treg$coefficients[parm] + lrint * z * se[parm],
                           extendInt = extendInt_upper, ...))
      return(c(lower$root, upper$root))
    }
    lims <- sapply(parm, limits)
    lower <- lims[1, ]
    upper <- lims[2, ]
  }

  # format output into data frame
  ci <- array(c(lower, upper), dim = c(length(parm), 2),
              dimnames = list(parm, pct))
  return(as.data.frame(ci))
}

#' @export
logLik.transreg <- function(treg) {
  # add degrees of freedom to allow AIC calculation
  logLik <- structure(treg$loglik, df = treg$df, class = "logLik")
  return(logLik)
}

#' @export
print.transreg <- function(treg) {
  cat("Call:\n")
  print(treg$call)
  cat("\n", "Coefficient estimates:\n")
  print(treg$coefficients)
  if (!is.null(treg$fixed)) {
    cat("\n", "Fixed coefficients:\n")
    print(treg$fixed)
  }
  if (treg$df == 1) {
    deg <- "degree"
  } else {
    deg <- "degrees"
  }
  cat("\n", "Log likelihood =", treg$loglik, "on", treg$df, deg,
      "of freedom.\n")
}

#' p-values for estimated parameters
#'
#' Calculates p-values for estimated parameters from a \code{transreg} model
#' using a normal approximation or a likelihood ratio chi-squared statistic.
#'
#' @param treg An object of class \code{transreg}.
#' @param parm A parameter name or vector of parameter names. If missing,
#'  p-values are calculated for all estimated parameters.
#' @param type The type of p-value. Options are \code{wald} for Wald p-values
#'  and \code{lr} for likelihood ratio p-values. The latter are more accurate
#'  but more computationally intensive.
#'
#' @return A named vector of p-values.
#'
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
pval.transreg <- function(treg, parm, type="wald") {
  # validate parameters
  if (missing(parm)) {
    parm <- names(treg$coefficients)
  } else if (anyNA(match(parm, names(treg$coefficients)))) {
    stop("Each name must match a parameter name.")
  }

  # determine type
  type_table <- c("wald", "lr")
  type <- type_table[pmatch(type, type_table)]
  if (is.na(type)) stop("Type not recognized.")

  # p-values
  if (type == "wald") {
    se <- sqrt(diag(treg$var)[parm])
    z <- abs(treg$coefficients[parm]) / se
    pvals <- 2 * pnorm(-z)
  } else {
    pval <- function(parm) {
      index <- match(parm, names(treg$coefficients))
      pvec <- treg$coefficients[-index]
      fixed <- c(treg$fixed, 0)
      names(fixed) <- c(names(treg$fixed), parm)
      parm_null <- stats::optim(pvec, treg$nlnL, fvec = fixed,
                                method = treg$optim_method)
      lnL_null <- -parm_null$value
      return(1 - pchisq(2 * (treg$loglik - lnL_null), 1))
    }
    pvals <- sapply(parm, pval)
  }
  return(pvals)
}

#' Summary of fitted transreg model
#'
#' Produces a summary of a fitted \code{transreg} model. The confidence level
#' and method for p-values and confidence intervals can be specified.
#'
#' @param treg An object of class \code{transreg}.
#' @param conf.level The confidence level (1 - \eqn{alpha}).
#' @param conf.type The type of confidence intervals and p-values. Options are
#'  \code{wald} for Wald and \code{lr} for likelihood ratio. This argument is
#'  passed to the \code{confint} and \code{pval} methods.
#'
#' @return A list with class \code{transreg_summary} that contains the
#'  following objects:
#'  \describe{
#'    \item{\code{call}}{The call to \code{transreg} with complete formal
#'      arguments.}
#'    \item{\code{dist}}{String naming the internal contact interval
#'      distribution.}
#'    \item{\code{lrt}}{A list containing the results of the global
#'      likelihood ratio test: \code{D} is the deviance, \code{df} is the
#'      degrees of freedom, \code{loglik} is the maximum log likelihood,
#'      \code{loglik_null} is the null log likelihood, and \code{p} is the
#'      p-value.}
#'    \item{\code{table}}{The coefficient table. Each row corresponds to an
#'      estimated parameter. The first column has point estimates, the second
#'      and third columns have confidence limits, and the last column has
#'      p-values.}
#'    \item{\code{type_name}}{String giving the method used to calculate
#'      p-values and confidence intervals.}
#'    \item{\code{xdist_name}}{A string giving the name of the external contact
#'      interval distribution; \code{NULL} if model has no \code{ext} term.}
#'  }
#'
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
summary.transreg <- function(treg, conf.level=0.95, conf.type="wald") {
  # get distribution names
  dist_names <- c("Exponential", "Log-logistic", "Weibull")
  index <- match(treg$dist, c("exponential", "loglogistic", "weibull"))
  dist_name <- dist_names[index]
  if (!is.null(treg$xdist)) {
    xdist_names <- c("Exponential", "Log-logistic", "Weibull")
    xindex <- match(treg$xdist, c("exponential", "loglogistic", "weibull"))
    xdist_name <- xdist_names[xindex]
  } else {
    xdist_name <- NULL
  }

  # fit with only fixed values
  if (treg$dist == "exponential") {
    pvec_null <- c("intercept" = 0)
  } else {
    pvec_null <- c("logshape" = 0, "intercept" = 0)
  }
  if (!is.null(treg$xdist)) {
    if (treg$xdist == "exponential") {
      pvec_null <- c("xintercept" = 0, pvec_null)
    } else {
      pvec_null <- c("xlogshape" = 0, "xintercept" = 0, pvec_null)
    }
  }
  pvec_null <- pvec_null[setdiff(names(pvec_null), names(treg$fixed))]
  df_null <- length(pvec_null)
  fit_null <- stats::optim(pvec_null, treg$nlnL, fvec = treg$fixed,
                           method = treg$optim_method)
  loglik_null <- -fit_null$value

  # get p-values and confidence limits
  if (!is.null(treg$coefficients)) {
    index <- match(conf.type, c("wald", "lr"))
    if (is.na(index)) stop("Confidence interval type not recognized.")
    type_name <- c("Wald", "Likelihood ratio")[index]
    table <- cbind(coef = treg$coefficients,
                   confint(treg, level = conf.level, type = conf.type),
                   p = pval(treg, type = conf.type))
  } else {
    table <- NULL
    type_name <- NULL
  }

  if (treg$df > df_null) {
    D <- 2 * (treg$loglik - loglik_null)
    df <- treg$df - df_null
    p <- 1 - pchisq(D, df)
    lrt <- list(D = D, df = df, loglik_null = loglik_null, p = p)
  } else if (treg$df > 0) {
    lrt <- "Null model"
  } else {
    lrt <- "All parameters fixed"
  }

  treg_summary <- structure(list(call = treg$call,
                                 dist_name = dist_name,
                                 fixed = treg$fixed,
                                 loglik = treg$loglik,
                                 lrt = lrt,
                                 table = table,
                                 type_name = type_name,
                                 xdist_name = xdist_name),
                            class = "transreg_summary")
  return(treg_summary)
}

#' @export
vcov.transreg <- function(treg) {
  return(treg$var)
}

# Methods for summary_transreg objects ====================================
# Need print method to prevent double printing without assignment.

#' Print summary of fitted transreg model
#'
#' Prints a summary of a fitted \code{transreg} model. The number of digits
#' to print for parameter estimates and p-values can be specified. The
#' p-values are formatted using \code{\link[base]{format.pval}}.
#'
#' @param treg_sum An object of class \code{transreg_summary}.
#' @param cdigits The minimum number of significant digits to print for
#'  parameter point and interval estimates.
#' @param pdigits The minimum number of significant digits to print for p-
#'  values. This is passed to \code{\link[base]{format.pval}}.
#'
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
print.transreg_summary <- function(treg_sum, cdigits=4, pdigits=3) {
  # print call, coefficients, p-values, and confidence limits
  cat("Call:\n")
  print(treg_sum$call)
  cat("\n")

  if (!is.null(treg_sum$table)) {
    cat("Contact intervals:", treg_sum$dist_name)
    if (!is.null(treg_sum$xdist)) {
      cat("\t(internal)\n", "                   ", treg_sum$xdist_name,
          "\t(external)\n", sep = "")
    }
    cat("\nConfidence intervals and p-values:", treg_sum$type_name, "\n")
    print(cbind(format(treg_sum$table[, 1:3], digits = cdigits),
                p = format.pval(treg_sum$table[, 4, drop = FALSE],
                                digits = pdigits)))
    cat("\n")
  }

  if (!is.null(treg_sum$fixed)) {
    cat("Fixed coefficients:\n")
    print(treg_sum$fixed)
    cat("\n")
  }

  # likelihood ratio test
  if (class(treg_sum$lrt) == "list") {
    lrt <- treg_sum$lrt
    cat("logLik(model) =", treg_sum$loglik, "\n")
    cat("logLik(null)  =", lrt$loglik_null, "\n")
    cat("  Chi^2 =", format(lrt$D, digits = cdigits), "on",
        lrt$df, "df:", "p =", format(lrt$p, digits = pdigits), "\n")
  } else {
    cat("logLik(model) =", treg_sum$loglik, "\n")
    cat(paste(" ", treg_sum$lrt, "\n"))
  }
}
