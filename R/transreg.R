#' Parametric regression models for infectious disease transmission
#'
#' Fits accelerated failure time models for infectious disease transmission 
#' data using right-censored data on contact intervals in ordered pairs of 
#' individuals. Who-infected-whom can be observed, unobserved, or partially 
#' observed.
#'
#' @param formula A formula of the form "response ~ terms". The response 
#'  must be an object returned by \code{\link[survival]{Surv}}. Only
#'  right-censored data is supported.
#' @param data A data frame containing the variables named in \code{formula}.
#' @param subset An expression indicating which rows of \code{data} should be
#'  included in the model fit.
#' @param na.action A missing-data filter applied to \code{data} after
#'  \code{subset}. Defaults to \code{options()$na.action}.
#' @param sus A string giving the name of the variable in \code{data} that
#'  contains the susceptible member of each pair.
#' @param dist A string partially matching a survival time distribution from 
#'  \code{\link{transreg.distributions}}. Current options are 
#'  \code{exponential}, \code{weibull} (the default), and \code{loglogistic}.
#' @param init A named vector of initial values for estimated coefficients.
#' @param fixed A named vector of fixed parameter values. These can include 
#'  terms in \code{formula}, \code{(Intercept)}, and \code{log(shape)}.
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
#'  \strong{Accelerated failure time models} are log-linear models for the 
#'  rate parameter \eqn{\lambda} in the specified survival time distribution. 
#'  Each coefficient can be interpreted as the log rate ratio for a one-unit 
#'  increase in the corresponding covariate, so each covariate has a 
#'  multiplicative effect on the rate. A single shape parameter \eqn{\gamma}
#'  is estimated for each stratum.
#' 
#' @return A list with class \code{transreg} that contains the following 
#'  objects:
#'  \describe{
#'    \item{\code{call}}{The call to \code{transreg} with complete formal 
#'      arguments.}
#'    \item{\code{coef}}{Named vector of estimated parameters.}
#'    \item{\code{df}}{The number of estimated coefficients.}
#'    \item{\code{dist}}{String naming the survival time distribution.}
#'    \item{\code{fixed}}{Named vector of fixed parameters.}
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
#'  }
#' 
#' @author Eben Kenah \email{ekenah@ufl.edu}
#' @references E Kenah (2011). Contact intervals, survival analysis of 
#'  epidemic data, and estimation of R_0. \emph{Biostatistics} 12(3): 
#'  548-566.
#' @references E Kenah and Y Sharker (2016)
#' @export
transreg <- function(
  formula, data, subset=NULL, na.action, sus, dist="weibull",
  init=NULL, fixed=NULL, optim_method="BFGS", ...
) {

  # fit accelerated failure time model using pairwise data

  # match arguments and ensure that formula argument is provided
  mcall <- match.call(expand.dots = FALSE)
  indx <- match(
    c("formula", "data", "subset", "na.action"), 
    names(mcall), nomatch = 0
  )
  if (indx[1] == 0) stop("A formula argument is required.")

  # find distribution from partial match
  dist_table <- c("exponential", "loglogistic", "lognormal", "weibull")
  dist <- dist_table[pmatch(dist, dist_table)]
  if (is.na(dist)) stop("Failure time distribution not recognized.")

  # pass model frame arguments to stats::model.frame
  model <- mcall[c(1, indx)]
  model[[1]] <- quote(stats::model.frame)
  mframe <- eval.parent(model)

  # get model responses and model matrix
  mterms <- attr(mframe, "terms")
  x <- model.matrix(mterms, mframe) 
  y <- model.response(mframe, "numeric")
  ymat <- data.frame(as.matrix(y))

  # factor of susceptibles in pairs with possible transmission
  if (missing(sus)) stop("Susceptible identifier not specified.")
  sus <- as.factor(data[ymat$status == 1, sus])

  # initial coefficient vector with log(shape) parameter
  beta <- rep(0, ncol(x))
  names(beta) <- colnames(x)
  total.infections <- nlevels(sus)
  total.time <- sum(ymat$time)
  beta["(Intercept)"] <- log(total.infections / total.time)
  if (dist == "exponential") {
    pvec <- beta
  } else {
    pvec <- c(0, beta)
    names(pvec) <- c("log(shape)", names(beta))
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

    # separate coefficients to be estimated and fixed coefficients
    fixed_indices <- match(names(fixed), names(pvec), nomatch = 0)
    if (any(fixed_indices == 0)) {
      stop("Each fixed value name must match a coefficient name.")
    }
    pvec <- pvec[-fixed_indices]
  }

  # full model negative log likelihood
  nlnL <- function(pvec, fvec=fixed) {
    pvec_lsindex <- match("log(shape)", names(pvec), nomatch = 0)
    fixed_lsindex <- match("log(shape)", names(fvec), nomatch = 0)
    if (pvec_lsindex > 0) {
        beta <- c(pvec[-pvec_lsindex], fvec)
        ymat$shape <- exp(pvec["log(shape)"])
    } else if (fixed_lsindex > 0) {
        beta <- c(pvec, fvec[-fixed_lsindex])
        ymat$shape <- exp(fvec["log(shape)"])
    } else {
        beta <- c(pvec, fvec)
    }
    ymat$rate <- exp(x[, names(beta), drop = FALSE] %*% beta)
    hmat <- subset(ymat, status == 1)
    return(transreg.nlnL(dist, hmat, ymat, sus))
  }

  # full model fit with point estimates and variance matrix
  fit <- stats::optim(
    pvec, nlnL, fvec = fixed, method = optim_method, hessian = TRUE, ...
  )
  if (length(fit$hessian) > 0) {
    coef <- fit$par
    var <- solve(fit$hessian)
  } else {
    coef <- NULL
    var <- NULL
  }

  output <- list(
    call = mcall,
    coef = coef,
    df = length(coef),
    dist = dist,
    fixed = fixed,
    loglik = -fit$value,
    model_matrix = x,
    nlnL = nlnL,
    optim_method = optim_method,
    response = y,
    sus = sus,
    var = var
  )
  class(output) <- "transreg"
  return(output)
}

# Internal methods used by transreg ---------------------------------------
# listed in alphabetical order

transreg.distributions <- list(
  exponential = list(
    haz = function(t, rate, shape) rate,
    cumhaz = function(t, rate, shape) rate * t
  ),
  loglogistic = list(
    haz = function(t, rate, shape) {
      shape * rate^shape * t^(shape - 1) / (1 + (rate * t)^shape)
    },
    cumhaz = function(t, rate, shape) log(1 + (rate * t)^shape)
  ),
  weibull = list(
    haz = function(t, rate, shape) shape * rate^shape * t^(shape - 1),
    cumhaz = function(t, rate, shape) (rate * t)^shape
  )
)

transreg.nlnL <- function(dist, hmat, ymat, sus) {
  haz <- transreg.distributions[[dist]]$haz
  cumhaz <- transreg.distributions[[dist]]$cumhaz

  # add hazards from all sources of infection for each infected susceptible
  hazards <- haz(t = hmat$time, rate = hmat$rate, shape = hmat$shape)
  hsums <- tapply(hazards, sus, sum)
  logh <- sum(log(hsums))

  # add cumulative hazards to get -log(survival)
  cumhazards <- cumhaz(t = ymat$time, rate = ymat$rate, shape = ymat$shape) 
  logS <- -sum(cumhazards)
  return(-logh - logS)
}

# Methods for transreg objects ============================================
# listed in alphabetical order

#' @export
coef.transreg <- function(treg) {
  return(treg$coef)
}

#' Confidence intervals for estimated parameters
#' 
#' Calculates confidence intervals for estimated parameters from a 
#' \code{transreg} model by inverting the Wald or likelihood ratio tests.
#' 
#' @param treg An object of class \code{treg}.
#' @param parm A parameter or vector of parameters. If missing, confidence 
#'  intervals are calculated for all estimated parameters.
#' @param level The confidence level (1 - \eqn{alpha}).
#' @param method The type of confidence interval. Current options are 
#'  \code{wald} for Wald confidence limits and \code{lr} for likelihood ratio
#'  confidence limits.
#' 
#' @return A data frame containing the lower and upper confidence limits with 
#'  column labeled with \eqn{\frac{\alpha}{2}} and 
#'  \eqn{1 - \frac{\alpha}{2}} expressed as percentages.
#' 
#' @author Eben Kenah \email{ekenah@ufl.edu}
#' @export
confint.transreg <- function(treg, parm, level=0.95, method="wald") {
  # validate parameters
  if (missing(parm)) { 
    parm <- names(treg$coef)
  } else if (anyNA(match(parm, names(treg$coef)))) {
    stop("Each parameter must match a coefficient name.")
  }

  # get percentages for lower and upper limits
  alpha <- 1 - level
  prb <- c(alpha / 2, 1 - alpha / 2)
  pct <- paste(signif(100 * prb, digits = 3), '%', sep = '') 

  # determine method
  method_table <- c("wald", "lr")
  method <- method_table[pmatch(method, method_table)]
  if (is.na(method)) stop("Method not recognized.")

  # find Z scores and standard error
  z <- qnorm(1 - alpha / 2)
  se <- sqrt(diag(treg$var)[parm])

  # confidence limits
  if (method == "wald") {
    lower <- treg$coef[parm] - z * se
    upper <- treg$coef[parm] + z * se
  } else {
    d <- qchisq(level, df = 1)
    limits <- function(parm) {
      pvec <- treg$coef[-match(parm, names(treg$coef))]
      parm_d <- function(val) {
        fixed <- c(treg$fixed, val)
        names(fixed) <- c(names(treg$fixed), parm)
        parm_fit <- stats::optim(
          pvec, treg$nlnL, fvec = fixed, method = treg$optim_method
        )
        return(2 * (treg$loglik + parm_fit$val) - d)
      }
      lower <- uniroot(parm_d, treg$coef[parm] + c(-1.5, -.5) * z * se[parm])
      upper <- uniroot(parm_d, treg$coef[parm] + c(.5, 1.5) * z * se[parm])
      return(c(lower$root, upper$root))
    }
    lims <- sapply(names(treg$coef), limits)
    lower <- lims[1, ]
    upper <- lims[2, ]
  }

  # format output into data frame
  ci <- array(
    c(lower, upper), dim = c(length(parm), 2), dimnames = list(parm, pct)
  )
  return(as.data.frame(ci))
}

#' @export
logLik.transreg <- function(treg) {
  # add degrees of freedom to allow AIC calculation
  logLik <- treg$loglik
  attr(logLik, 'df') <- treg$df
  return(logLik)
}

#' @export
print.transreg <- function(treg) {
  cat("Call:\n")
  print(treg$call)
  cat("\n", "Coefficient estimates:\n")
  print(treg$coef)
  if (!is.null(treg$fixed)) {
    cat("\n", "Fixed coefficients:\n")
    print(treg$fixed)
  }
  cat(
    "\n", "Log likelihood =", treg$loglik, "on", treg$df, 
    "degrees of freedom.\n"
  )
}

#' p-values for estimated parameters
#' 
#' Calculates p-values for estimated parameters from a \code{transreg} model 
#' using a normal approximation or a likelihood ratio chi-squared statistic.
#' 
#' @param treg An object of class \code{treg}.
#' @param parm A parameter or vector of parameters. If missing, p-values 
#'  are calculated for all estimated parameters.
#' @param method The type of p-value. Current options are \code{wald} for 
#'  normal approximation p-values and \code{lr} for likelihood ratio p-values.
#' 
#' @return A named vector of p-values.
#' 
#' @author Eben Kenah \email{ekenah@ufl.edu}
#' @export
pval.transreg <- function(treg, parm, method="wald") {
  # validate parameters
  if (missing(parm)) { 
    parm <- names(treg$coef)
  } else if (anyNA(match(parm, names(treg$coef)))) {
    stop("Each parameter must match a coefficient name.")
  }

  # determine method
  method_table <- c("wald", "lr")
  method <- method_table[pmatch(method, method_table)]
  if (is.na(method)) stop("Method not recognized.")

  # p-values
  if (method == "wald") {
    se <- sqrt(diag(treg$var)[parm])
    z <- abs(treg$coef[parm]) / se
    pvals <- 2 * pnorm(-z)
  } else {
    pval <- function(parm) {
      index <- match(parm, names(treg$coef))
      pvec <- treg$coef[-index]
      fixed <- c(treg$fixed, 0)
      names(fixed) <- c(names(treg$fixed), parm)
      parm_null <- stats::optim(
        pvec, treg$nlnL, fvec = fixed, method = treg$optim_method
      )
      lnL.null <- -parm_null$value
      return(1 - pchisq(2 * (treg$loglik - lnL.null), 1))
    }
    pvals <- sapply(parm, pval)
  }
  return(pvals)
}

#' Summary of fitted parametric regression model for ID transmission
#' 
#' Produces a summary of a fitted \code{transreg} model. The confidence level
#' and method for p-values and confidence intervals can be specified.
#' 
#' @param treg An object of class \code{transreg}.
#' @param conf.level The confidence level (1 - \eqn{alpha}).
#' @param conf.method The method for confidence intervals and p-values. 
#'  Current options are \code{wald} for Wald and \code{lr} for likelihood 
#'  ratio. These options correspond to those in the \code{confint} and 
#'  \code{pval} methods for \code{transreg} objects.
#' 
#' @return A list with class \code{summary_transreg} that contains the 
#'  following objects:
#'  \describe{
#'    \item{\code{call}}{The call to \code{transreg} with complete formal 
#'      arguments.}
#'    \item{\code{dist}}{String naming the survival time distribution.}
#'    \item{\code{lrt}}{A list containing the results of the global 
#'      likelihood ratio test: \code{D} is the deviance, \code{df} is the
#'      degrees of freedom, \code{loglik} is the maximum log likelihood, 
#'      \code{loglik_null} is the null log likelihood, and \code{p} is the 
#'      p-value.}
#'    \item{\code{table}}{The coefficient table. Each row corresponds to an 
#'      estimated parameter. The first column has point estimates, the second 
#'      and third columns have confidence limits, and the last column has 
#'      p-values.}
#'    \item{\code{type_name}}{A string giving the method used to calculate 
#'      p-values and confidence intervals.}
#'  }
#'
#' @author Eben Kenah \email{ekenah@ufl.edu}
#' @export
summary.transreg <- function(treg, conf.level=0.95, conf.type="wald") {
  # get pretty distribution and type name
  dist_names <- c("Exponential", "Log-logistic", "Lognormal", "Weibull")
  index <- match(
    treg$dist, 
    c("exponential", "loglogistic", "lognormal", "weibull")
  )
  dist_name <- dist_names[index]

  # fit with only fixed values
  if (treg$dist == "exponential") {
    pvec_null <- c("(Intercept)" = 0)
  } else {
    pvec_null <- c("log(shape)" = 0, "(Intercept)" = 0)
  }
  pvec_null <- pvec_null[setdiff(names(pvec_null), names(treg$fixed))]
  df_null <- length(pvec_null)
  fit_null <- stats::optim(
    pvec_null, treg$nlnL, fvec = treg$fixed, method = treg$optim_method
  )
  loglik_null <- -fit_null$value

  # make data  p-values and confidence limits
  index <- match(conf.type, c("wald", "lr"))
  if (is.na(index)) stop("Confidence interval type not recognized.")
  type_name <- c("Wald", "LR")[index]
  table <- cbind(
    coef = treg$coef, 
    confint(treg, level = conf.level, method = conf.type), 
    p = pval(treg)
  )

  if (treg$df > df_null) {
    D <- 2 * (treg$loglik - loglik_null)
    df <- treg$df - df_null
    p = 1 - pchisq(D, df)
    lrt <- list(D = D, df = df, p = p)
  } else {
    lrt <- NULL
  }

  treg_summary <- list(
    call = treg$call, 
    dist_name = dist_name, 
    fixed = treg$fixed,
    lrt = list(
      D = D, df = df, loglik = treg$loglik, loglik_null = loglik_null, p = p
    ),
    table = table,
    type_name = type_name
  )
  class(treg_summary) <- "summary_transreg"

  return(treg_summary)
}

#' @export
vcov.transreg <- function(treg) {
  return(treg$var)
}

# Methods for summary.transreg objects ====================================
# need print method to prevent double printing without assignment

#' Print summary of fitted parametric regression model for ID transmission
#' 
#' Prints a summary of a fitted \code{transreg} model. The number of digits
#' to print for parameter estimates and p-values can be specified. The 
#' p-values are formatted using \code{\link[base]{format.pval}}.
#' 
#' @param treg_sum An object of class \code{summary_transreg}.
#' @param cdigits The minimum number of significant digits to print for 
#'  parameter point and interval estimates.
#' @param pdigits The minimum number of significant digits to print for p-
#'  values. This is passed to \code{\link[base]{format.pval}}.
#' 
#' @author Eben Kenah \email{ekenah@ufl.edu}
#' @export
print.summary_transreg <- function(treg_sum, cdigits=4, pdigits=3) {
  # print call, coefficients, p-values, and confidence limits
  cat("Call:\n")
  print(treg_sum$call)
  cat("\n")

  cat(
    treg_sum$dist_name, " distribution estimates (", 
    treg_sum$type_name, "):\n", sep = ""
  )
  print(cbind(
    format(treg_sum$table[, 1:3], digits = cdigits),
    p = format.pval(treg_sum$table[, 4, drop = FALSE], digits = pdigits)
  ))
  cat("\n")

  if (!is.null(treg_sum$fixed)) {
    cat("Fixed coefficients:\n")
    print(treg_sum$fixed)
    cat("\n")
  }

  # likelihood ratio test
  lrt <- treg_sum$lrt
  cat("logLik(model)          =", lrt$loglik, "\n")
  cat("logLik(intercept only) =", lrt$loglik_null, "\n")
  cat(
    "  Chi^2 =", format(lrt$D, digits = cdigits), "on", 
    lrt$df, "df:", "p =", format(lrt$p, digits = pdigits), "\n"
  )
}

