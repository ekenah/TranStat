#' Nonparametric estimators for infectious disease transmission
#'
#' Calculates Nelson-Aalen or Kaplan-Meier curves for infectious disease 
#' transmission using right-censored and/or left truncated data on contact 
#' intervals in ordered pairs of individuals and infectious contacts from 
#' external sources with individuals.
#'
#' @param formula A formula of the form "response ~ terms". The response must 
#'  be an object returned by \code{\link[survival]{Surv}}. The terms are 
#'  converted into a factor with a level corresponding to each unique 
#'  combination of covariate values.
#' @param sus A string giving the name of the variable indicating the 
#'  susceptible member of each pair.
#' @param data A data frame containing the variables named in \code{formula}.
#' @param subset An expression indicating which rows of \code{data} should be
#'  included in the model fit. Default is to keep all rows.
#' @param na.action A missing-data filter applied to \code{data} after
#'  \code{subset}. Defaults to \code{options()$na.action}.
#' @param inithaz A named hazard function (up to a constant of proportionality) #'  to be used in the first iteration. The formal arguments in the function 
#'  need to match variable names \code{data} that appear in \code{formula}.
#' @param maxit The maximum number of EM algorithm iterations.
#' 
#' @author Eben Kenah \email{kenah.1@osu.edu}
#' @references E Kenah (2013). Nonparametric survival analysis of infectious 
#'  disease data. \emph{Journal of the Royal Statistical Society, Series B} 
#'  75(2): 277-303.

transfit <- function(
  formula, sus, data, subset=NULL, weights, na.action,
  inithaz=NULL, maxit=50
) {
  # calculate cumulative hazard or survival curves for transmission data

  # match arguments and ensure that formula argument is provided  
  mcall <- match.call(expand.dots = FALSE)
  indx <- match(
    c("formula", "data", "subset", "weights", "na.action"), 
    names(mcall), nomatch = 0
  )
  if (indx[1] == 0) stop("A formula argument is required.")

  # pass model frame arguments to stats::model.frame
  model <- mcall[c(1, indx)]
  model[[1]] <- quote(stats::model.frame)
  special <- "ext"
  model$formula <- terms(formula, special)
  mframe <- eval.parent(model)

  # get model matrix and responses
  mterms <- attr(mframe, "terms")
  if (any(attr(mterms, "order") > 1)) stop("Interactions not applicable.")
  attr(mterms, "intercept") <- 0  # supress intercept
  x <- model.matrix(mterms, mframe) 
  y <- model.response(mframe, "numeric")
  if (!survival::is.Surv(y)) stop("Response must be a Surv object.")
  ymat <- data.frame(as.matrix(y))
  attr(ymat, "type") <- attr(y, "type")

  # add indicator for external rows to response matrix
  if (is.null(attr(mterms, "specials")$ext)) {
    ext <- rep(0, nrow(x))
  } else {
    extname <- survival::untangle.specials(mterms, "ext")$vars
    ext <- x[, extname]
    if (sum(ext) == 0) {
      stop("Formula includes ext() term, but data has no external rows.")
    } 
  }
  ymat$ext <- ext

  # identify susceptibles in pairs with possible transmission
  if (missing(sus)) stop("Susceptible identifier not specified.")
  sus <- data[, sus]
  if (is.null(subset)) {
    ymat$sus <- sus
  } else {
    ymat$sus <- sus[eval(substitute(subset), data)]
  }

  # add initial hazard to response matrix
  if (is.null(inithaz)) {
    ymat$hazard <- ifelse(ymat$status, 1, 0)
  } else {
    ymat$hazard <- ifelse(ymat$status, eval(body(inithaz), data.frame(x)), 0)
  }

  # normalize hazard estimates
  ymat$pinfector <- 0
  for (sus in unique(ymat$sus[ymat$status == 1])) {
    hrows <- (ymat$sus == sus & ymat$status == 1)
    hazards <- ymat$hazard[hrows]
    hsum <- sum(hazards)
    ymat$pinfector[hrows] <- hazards / hsum
  }

  # convert model terms to factor
  if (ncol(x) == 0) {
    # ~1 on the right-hand side of the formula
    x <- factor(rep(1, nrow(x)))  
  } else {
    x <- survival::strata(as.data.frame(x))
  }


  # split responses into list with one element for each stratum
  ylist <- split(ymat, x)
#
#  mNA <- function(ystrat) {
#    # calculate number at risk as a function of analysis time
#    if (attr(ystrat, "type") == "right") {
#      stops <- sort(ystrat$time)
#      starts <- rep(0, nrow(ystrat))
#    } else if (attr(ystrat, "type") == "counting") {
#      stops <- sort(ystrat$stop)
#      starts <- sort(ystrat$start)
#    } else stop("Survival type not recognized.")
#    cumexits <- stepfun(stops, cumsum(c(0, table(stops))), right = TRUE)
#    cumentries <- stepfun(starts, cumsum(c(0, table(starts))), right = TRUE)
#
#    ctimes <- ystrat$stop[ystrat$status == 1]
#    nrisk <- cumexits(ctimes) - cumexits(ctimes)
#    Edeltas <- with(ystrat$pinfector[ystrat$status == 1])
#
#    # marginal Nelson-Aalen estimate
#    mNAincrements <- Edeltas / nrisk
#    mNAest <- stepfun(ctimes, cumsum(c(0, mNAincrements)))
#
#    # variance
#    mNAvar <- cumsum(Edeltas / nrisk^2)
#    mNAcov <- 0
#    susc <- ystrat$sus[ystrat$status == 1]
#    for (j in unique(susc)) {
#      mNAcovj <- cumsum(mNAincrements * (susc == j))
#      mNAcov <- mNAcov + mNAcovj^2
#    }
#    mNAvar <- stepfun(ctimes, cumsum(c(0, 2 * mNAvar - mNAcov)))
#
#    return(list(est = mNAest, var = mNAvar))
#  }
#
#  mNAlist <- lapply(ylist, mNA)
#
#
#
#
#  # more iterations if WIW not completely observed and maxit > 1
#
#  if (length(mNAlist) == 1) mNA = mNAlist[[1]]

  return(list(
    x = x,
    y = y
  ))
  
}
