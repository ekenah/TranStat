#' Semiparametric regression models for infectious disease transmission
#'
#' Fits semiparametric pairwise regression models for infectious disease   
#' transmission using right-censored and/or left-truncated data on contact 
#' intervals in ordered pairs of individuals and infectious contact from 
#' external sources with individuals. Uses the Cox relative risk function. 
#' External source of infection are handled by stratifying on an external 
#' row indicator.
#'
#' @param formula A formula of the form "response ~ terms". The response 
#'  must be an object returned by \code{\link[survival]{Surv}} of type 
#'  "right" or "counting". The formula can include \code{offset} terms for 
#'  fixed coefficients, and an external row indicator can be included in a 
#'  \code{strata} term. 
#' @param sus A string giving the name of the variable in \code{data} that
#'  contains the susceptible member of each pair.
#' @param data A data frame containing the variables named in \code{formula}.
#' @param weights A vector of infector probabilities for the data rows. If 
#'  who-infected-whom is not completely observed, these will be updated in an
#'  EM algorithm. If missing, all possible infectors of each susceptible are
#'  given equal initial weights.
#' @param subset An expression indicating which rows of \code{data} should be
#'  included in the model fit.
#' @param na.action A missing-data filter applied to \code{data} after
#'  \code{subset}. Defaults to \code{options()$na.action}.
#' @param itermax The maximum number of EM algorithm iterations if who-
#'  infected-whom is not observed.
#' @param degf The degrees of freedom for the smoothing spline used to obtain 
#'  the hazard function(s) from the estimated cumulative hazard(s). The default 
#'  is to use twice the natural logarithm of the total number of transmissions.
#' @param L1tol The EM algorithm stops when the L1 distance between the old and 
#'  new weights is below \code{L1tol} per transmission event.
#' @param ... Further arguments to be passed to \code{\link[survival]{coxph}},
#'  such as \code{init}, \code{ties}, or \code{control}.
#' 
#' @return A list of class \code{transph} with the following components:
#'  \describe{
#'    \item{\code{call}}{The call to \code{transph} with complete formal 
#'      arguments.}
#'    \item{\code{coefficients}}{The estimated coefficients.}
#'    \item{\code{coxph}}{The \code{\link[survival]{coxph}} object returned by 
#'      the final Cox regression.}
#'    \item{\code{df}}{The number of estimated coefficients.}
#'    \item{\code{formula}}{The model formula.}
#'    \item{\code{iter}}{The number of iterations. It is one for a model in 
#'      which who-infects-whom is completely observed. Otherwise, it is the 
#'      number of EM algorithm iterations (including the initial model).}
#'    \item{\code{L1tol}}{The L1 tolerance per transmission event used to halt 
#'      the EM algorithm.}
#'    \item{\code{loglik}}{The \code{loglik} element in the 
#'      \code{\link[survival]{coxph.object}} returned by the final Cox 
#'      regression, which contains the log partial likelihood at the 
#'      initial coefficient values and the fitted coefficient values. For a 
#'      null model, only the first element is present.}
#'    \item{\code{method}}{The \code{method} element in the 
#'      \code{\link[survival]{coxph.object}} returned by the final Cox 
#'      regression, which is the name of approximation used to handle ties.}
#'    \item{\code{spline_df}}{The degrees of freedom in the smoothing spline 
#'      used to calculate the contact interval hazard function(s).}
#'    \item{\code{sus}}{The vector identifying the susceptible member of each #'      pair in the data used to fit the model.}
#'    \item{\code{var}}{The covariance matrix for the coefficient estimates. #'      The rows and columns are named for the corresponding covariates.}
#'    \item{\code{weights}}{The final vector of estimated infector 
#'      probabilities.}
#'  }
#' 
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
transph <- function(formula, sus, data, weights, subset, na.action, degf,
                    itermax=25, L1tol=1e-04, ...)
{
  # fit semiparametric model using pairwise data
  
  # match arguments and ensure that formula argument is provided
  mcall <- match.call(expand.dots = FALSE)
  indx <- match(
    c("formula", "data", "subset", "na.action"), 
    names(mcall), nomatch = 0
  )
  if (indx[1] == 0) stop("A formula argument is required.")

  # pass model frame arguments to stats::model.frame
  model <- mcall[c(1, indx)]
  model[[1]] <- quote(stats::model.frame)
  special <- c("strata")
  model$formula <- terms(formula, special)
  mframe <- eval.parent(model)

  # get model matrix and responses
  mterms <- attr(mframe, "terms")
  xmat <- data.frame(model.matrix(mterms, mframe))
  y <- model.response(mframe, "numeric")
  ymat <- data.frame(as.matrix(y))
  attr(ymat, "type") <- attr(y, "type")

  # identify susceptibles in pairs with possible transmission
  if (missing(sus)) stop("Susceptible identifier not specified.")
  if (class(substitute(sus)) == "character") sus <- as.name(sus)
  if (missing(subset)) {
    subset <- NULL
    sus <- eval(substitute(sus), data)
  } else {
    sus <- sus[eval(substitute(subset), data)]
  }

  # determine whether who-infected-whom is observed
  Vjsize <- tapply(ymat$status, sus, sum)   # ymat$status is always 0/1
  Vjmax <- max(Vjsize)
  Vjsize_sus <- sapply(sus, function(j) Vjsize[as.character(j)])
  ntrans <- length(unique(sus[ymat$status == 1]))

  # prepare data for weighted Cox regression
  if (Vjmax > 1) {
    # who-infected-whom not completely observed
    if (missing(degf)) {
      # determine default degrees of freedom for smoothing spline
      degf <- ceiling(2 * log(ntrans))
    } 

    # expand data to include copies of rows with possible transmission
    copyrows <- (ymat$status == 1 & Vjsize_sus > 1)
    data0 <- data[copyrows,]
    cdata <- rbind(data, data0)

    # expand y to include censored copies of outcomes
    ymat0 <- ymat[copyrows,]
    ymat0$status <- 0

    # rebuild y as a Surv object cy
    cy <- as.matrix(rbind(ymat, ymat0))
    attr(cy, "dim") <- dim(cy)
    attr(cy, "row.names") <- NULL
    attr(cy, "type") <- attr(ymat, "type")
    dimnames(cy) <- list(c(rownames(ymat), rownames(ymat0)), 
                         names(ymat))
    attr(cy, "class") <- "Surv"
  } else {
    # who-infected-whom observed
    cdata <- data
    cy <- y
  }

  # weighted Cox regression and EM algorithm (if needed)
  environment(formula) <- environment() # correct for NSE of model.frame
  cformula <- update(formula, cy ~ .)

  creg <- NULL  # shield from any "creg" defined outside the function
  L1dist <- Inf
  iter <- 0
  while (iter < itermax) {
    if (iter > 0) {
      creg0 <- creg
    }

    if (Vjmax > 1) {
      if (missing(weights)) {
        # give all possible infectors equal weight
        weights <- rep(1, nrow(data))
        weights[ymat$status == 1] <- 1 / Vjsize_sus[ymat$status == 1]
      } else if (!is.null(creg)) {
        # re-weight possible infectors using previous Cox regression
        old_weights <- weights
        weights <- transph.weights(creg, data, ymat, sus, degf)

        # calculate mean L1 distance between old and new weights
        L1dist <- sum(abs(weights - old_weights)) / ntrans
      }

      # construct cweights for new Cox regression
      weights0 <- 1 - weights[copyrows]
      cweights <- c(weights, weights0)

    } else {
      cweights <- weights
    }

    creg <- survival::coxph(cformula, cdata, cweights, subset, 
                            na.action, ...)   
    if (Vjmax <= 1 | L1dist < L1tol) {
      break
    } else {
      iter <- iter + 1
    }
  }
  if (Vjmax > 1 & iter > itermax) {
    warning("Iteration limit reached without specified L1 tolerance.")
  }

  # covariance matrix for coefficient estimates
  var <- creg$var
  if (Vjmax > 1 & length(coef(creg)) > 0) {
    # who-infected-whom not completely observed and non-null model
    csus <- c(sus, sus[copyrows])
    cxmat <- rbind(xmat, xmat[copyrows,])
    cymat <- data.frame(as.matrix(cy), row.names = NULL)
    attr(cymat, "type") <- attr(cy, "type")

    var <- solve(transph.information(creg, cdata, cymat, cxmat, csus))
  }
  rownames(var) <- names(creg$coefficients)
  colnames(var) <- names(creg$coefficients)

  # define transph object
  output <- structure(list(call = mcall,
                           coefficients = creg$coefficients,
                           coxph = creg,
                           df = length(creg$coefficients),
                           formula = formula,
                           iter = iter,
                           L1tol = L1tol,
                           loglik = creg$loglik,
                           method = creg$method,
                           nevent = ntrans,
                           spline_df = degf,
                           sus = sus,
                           var = var,
                           weights = weights),
                      class = "transph")
  return(output)
}

# Internal methods used by transph -------------------------------------------
# listed in alphabetical order

transph.weights <- function(creg, data, ymat, sus, degf) 
{
  # calculate baseline hazard estimates
  mbreslow <- survival::survfit(creg)

  # get times for possible transmission
  hmat <- subset(ymat, status == 1)
  if (attr(ymat, "type") == "right") {
    times <- hmat$time
  } else if (attr(ymat, "type") == "counting") {
    times <- hmat$stop
  } else {
    stop(paste("Unsupported Surv type: ", attr(ymat, "type")))
  }
  events <- ymat$status == 1

  # smoothed "baseline" cumulative hazard estimates
  if ("strata" %in% names(mbreslow)) {
    basehazards <- rep(0, nrow(hmat))

    # test for stratum membership
    stest <- gsub(",", "&", gsub("=", "==", names(mbreslow$strata)))

    # iterate through strata
    end <- 0
    for (s in 1:length(mbreslow$strata)) {
      start <- end + 1
      end <- end + mbreslow$strata[s]

      # predict baseline hazards within stratum
      n.event <- mbreslow$n.event[start:end]
      mbtimes <- mbreslow$time[start:end]
      mbsurv <- mbreslow$surv[start:end]
      mbstderr <- mbreslow$std.err[start:end]
      smooth_cumhaz <- smooth.spline(mbtimes[n.event > 0], 
                                     -log(mbsurv[n.event > 0]),
                                     w = mbstderr[n.event > 0]^(-2), df = degf)

      # put baseline hazards into stratum elements of vector
      stratum <- with(data[events,], eval(parse(text = stest[s])))
      bhazards <- predict(smooth_cumhaz, times[stratum], deriv = 1)$y
      basehazards[stratum] <- bhazards
    }
  } else {
    mbtimes <- with(mbreslow, time[n.event > 0])
    mbchaz <- with(mbreslow, -log(surv[n.event > 0]))
    mbweights <- with(mbreslow, std.err[n.event > 0]^(-2))
    smooth_cumhaz <- smooth.spline(mbtimes, mbchaz, w = mbweights, df = degf)
    basehazards <- predict(smooth_cumhaz, times, deriv = 1)$y
  }

  # hazard ratios
  hdata <- data[events,]
  if (is.null(coef(creg))) {
    hazratios <- 1
  } else {
    hazratios <- predict(creg, hdata, type = "risk", reference = "sample")
  }

  # recalculate weights for possible infectors using new hazards
  hazards <- basehazards * hazratios
  hsus <- sus[events]
  hsums <- tapply(hazards, hsus, sum)
  hweights <- hazards / sapply(hsus, function(j) hsums[as.character(j)])

  # create complete vector of weights
  weights <- rep(1, nrow(data))
  weights[events] <- hweights

  return(weights)
}

transph.information <- function(creg, cdata, cymat, cxmat, csus) 
{
  # get information matrix from Cox regression
  Imat <- solve(vcov(creg))

  # get times of possible transmission
  if (attr(cymat, "type") == "right") {
    times <- cymat$time
  } else if (attr(cymat, "type") == "counting") {
    times <- cymat$stop
  } else {
    stop(paste("Unsupported Surv type: ", attr(cymat, "type")))
  }
  events <- cymat$status == 1

  cdetail <- coxph.detail(creg)
  if ("strata" %in% names(cdetail)) {
    # test for stratum membership
    stest <- gsub(",", "&", gsub("=", "==", names(cdetail$strata)))

    # iterate through strata
    Uij <- matrix(rep(0, sum(cymat$status == 1) * length(names(coef(creg)))),
                  ncol = length(names(coef(creg))))
    end <- 0
    for (s in 1:length(cdetail$strata)) {
      start <- end + 1
      end <- end + cdetail$strata[s]

      # calculate Schoenfeld-type residuals for each event pair in stratum
      stratum <- with(cdata, eval(parse(text = stest[s])))
      stimes <- cdetail$time[start:end]
      eindx <- match(times[stratum & events], stimes) + start - 1
      Uij_diff <- (as.matrix(cxmat[stratum & events, names(coef(creg))]) 
                   - as.matrix(cdetail$means)[eindx,])
      Uij_wts <- creg$weights[stratum & events]
      Uij[stratum[events],] <- Uij_diff * Uij_wts  # order important
    }
  } else {
    # calculate Schoenfeld-type residual for each pair ij with an event
    eindx <- match(times[events], cdetail$time)
    Uij_diff <- (cxmat[events, names(coef(creg))] 
                 - as.matrix(cdetail$means)[eindx,])
    Uij_wts <- creg$weights[seq(1, nrow(cymat))[events]]
    Uij <- as.matrix(Uij_diff * Uij_wts)   # order important due to recycling
  }

  # calculate sum of residuals for each susceptible j
  Uj <- do.call(rbind, by(Uij, csus[events], colSums, simplify = FALSE))

  Uij2 <- t(Uij) %*% Uij
  Uj2 <- t(Uj) %*% Uj

  return(Imat - Uij2 + Uij2) 
}  


# Methods for transph objects ================================================

# define transph.object, anova.transph, logLik.transph, predict.transph, 
# residuals.transph, summary.transph, survfit.transph

#' Confidence intervals for estimated coefficients
#' 
#' Calculates confidence intervals for coefficient estimates from a 
#' \code{transph} model by inverting the Wald or likelihood ratio tests.
#' 
#' @param treg An object of class \code{transph}.
#' @param parm A coefficient name or vector of coefficient names. If missing, #'  confidence intervals are calculated for all estimated coefficients.
#' @param level The confidence level (1 - \eqn{alpha}).
#' @param type The type of confidence interval. Current options are \code{wald} 
#'  for Wald confidence limits and \code{lr} for likelihood ratio confidence 
#'  limits. The latter are more accurate but more computationally intensive.
#' 
#' @return A data frame with one row for each coefficient in \code{parm} that 
#'  contains the lower and upper confidence limits in columns labeled 
#'  \eqn{\frac{\alpha}{2}} and \eqn{1 - \frac{\alpha}{2}} expressed as 
#'  percentiles.
#' 
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
confint.transph <- function(creg, parm, level=0.95, type="wald") 
{
  if (missing(parm)) { 
    parm <- names(creg$coefficients)
  } else if (anyNA(match(parm, names(creg$coefficients)))) {
    stop("Each parameter must match a coefficient name.")
  }

  # get percentiles for lower and upper limits
  alpha <- 1 - level
  prb <- c(alpha / 2, 1 - alpha / 2)
  pct <- paste(signif(100 * prb, digits = 3), '%', sep = '') 

  # determine confidence interval type
  type_table <- c("wald", "lr")
  type <- type_table[pmatch(type, type_table)]
  if (is.na(type)) stop("Type not recognized.")

  # find Z scores and standard error
  z <- qnorm(1 - alpha / 2)
  se <- sqrt(diag(creg$var)[parm])

  # confidence limits
  if (type == "wald") {
    # Wald confidence limits
    lower <- creg$coefficients[parm] - z * se
    upper <- creg$coefficients[parm] + z * se
  } else {
    # likelihood ratio confidence limits
    d <- qchisq(level, df = 1)
    limits <- function(parm) 
    {
      parm_d <- function(val) 
      {
        # generate character vector to update model formula
        parm_offset <- paste("~ . - ", parm, " + offset(",
                             as.character(val), " * ", parm, ")", sep = "")
        pformula <- update(creg$formula, parm_offset)

        # call transph with parm coefficient fixed and itermax = 1
        pargs <- as.list(creg$call)[2:length(creg$call)]
        pargs$formula <- pformula
        pargs$init <- creg$coefficients[-match(parm, names(creg$coefficients))]
        pargs$weight <- creg$weights
        pargs$itermax <- 1  # keep estimated weights from full model
        parm_fit <- do.call(transph, pargs)
        lnL_model <- creg$loglik[2]
        lnL_null <- parm_fit$loglik[length(parm_fit$loglik)]
        return(2 * (lnL_model - lnL_null) - d)
      }
      lower <- list(root = -Inf)
      lower <- uniroot(
        parm_d, creg$coefficients[parm] + c(-2, -.5) * z * se[parm],
        extendInt = "downX"
      )
      upper <- list(root = Inf)
      try(upper <- uniroot(
        parm_d, creg$coefficients[parm] + c(.5, 2) * z * se[parm],
        extendInt = "upX"
      ))
      return(c(lower$root, upper$root))
    }
    lims <- sapply(parm, limits)
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
logLik.transph <- function(creg) 
{
  indx <- length(creg$loglik)
  logLik <- structure(creg$loglik[index], df = creg$df, class = "logLik")
  return(logLik)
}

#' @export
print.transph <- function(creg) 
{
  cat("Call:\n")
  print(creg$call)
  cat("\n", "Coefficient estimates:\n")
  print(creg$coefficients)
  if (creg$df == 1) {
    deg <- "degree"
  } else {
    deg <- "degrees"
  }
  cat("\n", "Log likelihood =", creg$loglik, "on", creg$df, deg, 
      "of freedom.\n")  
}

#' p-values for estimated coefficients
#' 
#' Calculates p-values for coefficient estimates from a \code{transph} model 
#' using a normal approximation or a likelihood ratio chi-squared statistic.
#' 
#' @param treg An object of class \code{transph}.
#' @param parm A coefficient name or vector of coefficient names. If missing, 
#'  p-values are calculated for all estimated parameters.
#' @param type The type of p-value. Options are \code{wald} for Wald p-values 
#'  and \code{lr} for likelihood ratio p-values. The latter are more accurate 
#'  but more computationally intensive.
#' 
#' @return A named vector of p-values.
#' 
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
pval.transph <- function(creg, parm, type="wald") 
{
  # validate parameters
  if (missing(parm)) { 
    parm <- names(creg$coefficients)
  } else if (anyNA(match(parm, names(creg$coefficients)))) {
    stop("Each name must match a coefficient name.")
  }

  # determine type
  type_table <- c("wald", "lr")
  type <- type_table[pmatch(type, type_table)]
  if (is.na(type)) stop("Type not recognized.")

  # p-values
  if (type == "wald") {
    se <- sqrt(diag(creg$var)[parm])
    z <- abs(creg$coefficients[parm]) / se
    pvals <- 2 * pnorm(-z)
  } else {
    pval <- function(parm) 
    {
      # generate character vector to update model formula
      parm_null <- paste("~ . -", parm)
      pformula <- update(creg$formula, parm_null)

      # call transph without parm
      pargs <- as.list(creg$call)[2:length(creg$call)]
      pargs$formula <- pformula
      pargs$init <- creg$coefficients[-match(parm, names(creg$coefficients))]
      pargs$weight <- creg$weights
      pargs$itermax <- 1  # keep estimated weights from full model
      parm_fit <- do.call(transph, pargs)

      lnL_model <- creg$loglik[2]
      lnL_null <- parm_fit$loglik[length(parm_fit$loglik)]
      return(1 - pchisq(2 * (lnL_model - lnL_null), 1))
    }
    pvals <- sapply(parm, pval)
  }
  return(pvals)
}

#' Summary of fitted transph model
#' 
#' Produces a summary of a fitted \code{transreg} model. The confidence level
#' and method for p-values and confidence intervals can be specified.
#' 
#' @param creg An object of class \code{transph}.
#' @param conf.level The confidence level (1 - \eqn{alpha}).
#' @param conf.type The type of confidence intervals and p-values. Options are 
#'  \code{wald} for Wald and \code{lr} for likelihood ratio. This argument is 
#'  passed to the \code{confint} and \code{pval} methods.
#' 
#' @return A list of class \code{transph_summary} with the following 
#'  components:
#'  \describe{
#'    \item{\code{call}}{The call to \code{transph} with complete formal 
#'      arguments.}
#'    \item{\code{iter}}{The number of iterations used to fit \code{creg}.}
#'    \item{\code{loglik}}{The \code{loglik} element of the \code{transph} 
#'      object.}
#'    \item{\code{lrt}}{A list containing the results of the global 
#'      likelihood ratio test: \code{D} is the deviance, \code{df} is the
#'      degrees of freedom, \code{loglik} is the maximum log likelihood, 
#'      \code{loglik_null} is the null log likelihood, and \code{p} is the 
#'      p-value.}
#'    \item{\code{nevent}}{The number of transmissions in the data used to fit 
#'      \code{creg}.}
#'    \item{\code{table}}{The coefficient table. Each row corresponds to an 
#'      estimated parameter. The first column has point estimates, the second 
#'      and third columns have confidence limits, and the last column has 
#'      p-values.}
#'    \item{\code{type_name}}{String giving the method used to calculate 
#'      p-values and confidence intervals.}
#'  }
#'
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
summary.transph <- function(creg, conf.level=0.95, conf.type="wald", 
                            cdigits=4, pdigits=3) 
{
  # determine type
  type_table <- c("wald", "lr")
  conf.type <- type_table[pmatch(conf.type, type_table)]
  if (is.na(conf.type)) stop("Type not recognized.")  

  # likelihood ratio test
  if (!is.null(creg$coefficients)) {
    D <- 2 * with(creg, loglik[2] - loglik[1])
    df <- creg$df
    p <- 1 - pchisq(D, df)
    lrt <- list(D = D, df = df, p = p)
  } else {
    lrt <- "Null model"
  }

  # p-values and confidence limits
  if (!is.null(creg$coefficients)) {
    index <- match(conf.type, c("wald", "lr"))
    if (is.na(index)) stop("Confidence interval type not recognized.")
    type_name <- c("Wald", "Likelihood ratio")[index]
    table <- cbind(coef = creg$coefficients, 
                   confint(creg, level = conf.level, type = conf.type), 
                   p = pval(creg, type = conf.type))
  } else {
    table <- NULL
    type_name <- NULL
  }  

  creg_summary <- structure(list(call = creg$call, 
                                 iter = creg$iter,
                                 loglik = creg$loglik,
                                 lrt = lrt,
                                 nevent = creg$nevent,
                                 table = table,
                                 type_name = type_name),
                            class = "transph_summary")
  return(creg_summary)  
}

#' Compute a survival curve from a transph model
#'
#' Computes
survfit.transph <- function(treg, newdata, se.fit=TRUE, conf.int=0.95,
                            type, vartype, conf.type, censor=FALSE,
                            start.time, na.action) 
{
  # create survival or cumulative hazard curves
  if (class(treg) != "transph") {
    stop("First argument should be a transph object.")
  }
  if (missing(newdata)) stop("No newdata argument provided")
}

#' @export
vcov.transph <- function(creg) 
{
  return(creg$var)
}

# Methods for summary_transph objects =====================================
# Need print method to prevent double printing without assignment.

#' Print summary of fitted transph model
#' 
#' Prints a summary of a fitted \code{transph} model. The number of digits
#' to print for parameter estimates and p-values can be specified. The 
#' p-values are formatted using \code{\link[base]{format.pval}}.
#' 
#' @param treg_sum An object of class \code{transph_summary}.
#' @param cdigits The minimum number of significant digits to print for 
#'  parameter point and interval estimates.
#' @param pdigits The minimum number of significant digits to print for p-
#'  values. This is passed to \code{\link[base]{format.pval}}.
#' 
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @export
print.transph_summary <- function(creg_sum, cdigits=4, pdigits=3)
{
  # print call, iter, and nevent
  cat("Call:\n")
  print(creg_sum$call)
  cat("Transmission events: ", creg_sum$nevent, "\n")
  cat("EM iterations: ", creg_sum$iter, "\n")
  cat("\n")

  # print table of coefficient estimates, confidence limits, and p-values
  if (!is.null(creg_sum$table)) {
    cat("\nConfidence intervals and p-values:", creg_sum$type_name, "\n")
    print(cbind(format(creg_sum$table[, 1:3], digits = cdigits), 
                p = format.pval(creg_sum$table[, 4, drop = FALSE], 
                                digits = pdigits)))
    cat("\n")
  }

  # print global likelihood ratio test
  if (class(creg_sum$lrt) == "list") {
    lrt <- creg_sum$lrt
    cat("logLik(model) =", creg_sum$loglik[2], "\n")
    cat("logLik(null)  =", creg_sum$loglik[1], "\n")
    cat("  Chi^2 =", format(lrt$D, digits = cdigits), "on", 
        lrt$df, "df:", "p =", format(lrt$p, digits = pdigits), "\n")
  } else {
    cat("logLik(model) =", creg_sum$loglik[1], "\n")
    cat(paste(" ", creg_sum$lrt, "\n"))
  }  
}
