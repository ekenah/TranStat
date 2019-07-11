#' Semiparametric regression models for infectious disease transmission
#'
#' Fits semiparametric pairwise regression models for infectious disease   
#' transmission using right-censored and/or left-truncated data on contact 
#' intervals in ordered pairs of individuals and infectious contact from 
#' external sources with individuals. Uses the Cox relative risk function.
#'
#' @param formula A formula of the form "response ~ terms". The response 
#'  must be an object returned by \code{\link[survival]{Surv}} of type "right" 
#'  or "counting". The formula can include \code{offset} terms for fixed 
#'  coefficients, and the external row indicator must be included in a 
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
#'  such as \code{init}, \code{ties}, and \code{control}.


transph <- function(formula, sus, data, weights=NULL, subset=NULL, na.action, 
                    itermax=25, degf=NULL, L1tol=1e-05, ...)
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
  x <- model.matrix(mterms, mframe)   
  y <- model.response(mframe, "numeric")
  ymat <- data.frame(as.matrix(y))
  attr(ymat, "type") <- attr(y, "type")

  # identify susceptibles in pairs with possible transmission
  if (missing(sus)) stop("Susceptible identifier not specified.")
  if (class(substitute(sus)) == "character") sus <- as.name(sus)
  if (missing(subset)) {
    sus <- eval(substitute(sus), data)
  } else {
    sus <- sus[eval(substitute(subset), data)]
  }

  # determine whether who-infected-whom is observed
  Vjsize <- tapply(ymat$status, sus, sum)   # ymat$status is always 0/1
  Vjmax <- max(Vjsize)
  Vjsize_sus <- sapply(sus, function(j) Vjsize[as.character(j)])
  ntrans <- length(unique(sus[ymat$status == 1]))

  if (Vjmax > 1) {
    # who-infected-whom is not completely observed
    if (is.null(degf)) {
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
    # who-infected-whom is observed
    cdata <- data
    cy <- y
  }

  # weighted Cox regression
  environment(formula) <- environment() # correct for NSE of model.frame
  cformula <- update(formula, cy ~ .)

  L1dist <- Inf
  iter <- 0
  while (iter <= itermax) {
    if (iter > 0) {
      creg0 <- creg
    }

    if (Vjmax > 1) {
      if (is.null(weights)) {
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

  # calculate covariance matrix
  beta <- coef(creg)
  if (!is.null(beta)) {
    creg_details <- survival::coxph.detail(creg)
  }

  # return transph object
  return(creg)
}

# Internal methods used by transph -------------------------------------------
# listed in alphabetical order

transph.weights <- function(creg, data, ymat, sus, degf) 
{
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

  # smoothed "baseline" cumulative hazard estimates
  if ("strata" %in% names(mbreslow)) {
    basehazards <- rep(0, nrow(hmat))

    # test for stratum membership
    stest <- gsub(",", "&", gsub("=", "==", names(mbreslow$strata)))

    # iterate through strata
    start <- 0
    end <- 0
    for (s in 1:length(mbreslow$strata)) {
      end <- start + mbreslow$strata[s]
      n.event <- mbreslow$n.event[start:end]
      mbtimes <- mbreslow$time[start:end]
      mbchaz <- -log(mbreslow$surv[start:end])
      mbstderr <- mbreslow$std.err[start:end]
      smooth_cumhaz <- smooth.spline(mbtimes[n.event > 0], mbchaz[n.event > 0],
                                     w = mbstderr[n.event > 0]^(-2))
      bhazards <- predict(smooth_cumhaz, times, deriv = 1)$y
      stratum <- with(data, parse(text = stest[s]))
      basehazards[stratum] <- bhazards
      start <- end + 1
    }
  } else {
    mbtimes <- with(mbreslow, time[n.event > 0])
    mbchaz <- with(mbreslow, -log(surv[n.event > 0]))
    mbweights <- with(mbreslow, std.err[n.event > 0]^(-2))
    smooth_cumhaz <- smooth.spline(mbtimes, mbchaz, w = mbweights, df = degf)
    basehazards <- predict(smooth_cumhaz, times, deriv = 1)$y
  }

  # hazard ratios
  hdata <- data[ymat$status == 1,]
  if (is.null(coef(creg))) {
    hazratios <- 1
  } else {
    hazratios <- predict(creg, hdata, type = "risk", reference = "sample")
  }

  # recalculate weights for possible infectors using new hazards
  hazards <- basehazards * hazratios
  hsus <- sus[ymat$status == 1]
  hsums <- tapply(hazards, hsus, sum)
  hweights <- hazards / sapply(hsus, function(j) hsums[as.character(j)])

  # create complete vector of weights
  weights <- rep(1, nrow(data))
  weights[ymat$status == 1] <- hweights

  return(weights)
}

transph.covariance <- function() 
{
}


# Methods for transph objects ================================================

# define transph.object, anova.transph, logLik.transph, predict.transph, 
# residuals.transph, summary.transph, survfit.transph

confint.transph <- function(creg, parm, level=0.95, type="wald") 
{
}

logLik.transph <- function(creg) 
{
}

print.transph <- function(creg) 
{
}

pval.transph <- function(creg, parm, type="wald") 
{
}

summary.transph <- function(creg, conf.level=0.95, conf.type="wald", 
                            cdigits=4, pdigits=3) 
{
}

survfit.transph <- function(formula, ...) 
{
  # create survival or cumulative hazard curves
}

vcov.transph <- function(creg) 
{
}

