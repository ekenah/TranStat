## roxygen tags requested by devtools::use_rcpp()
#' @useDynLib transtat
#' @importFrom Rcpp sourceCpp

## Make survival::Surv available to transtat users
# This is a hack; using importFrom(survival, Surv) didn't work.
#' @export
Surv <- survival::Surv

## Generic method for p-values (similar to coef or confint)
#' p-values for model parameters.
#' 
#' Calculates p-values for one or more estimated parameters.
#' 
#' @param object A fitted model object.
#' @param parm A parameter or vector of parameters. If missing, p-values are
#'  calculated for all estimated parameters.
#' @param ... Additional arguments for methods.
#'
#' @details
#'  \code{pval} is a generic function. The default method calculates Wald
#'  confidence intervals and assumes that \code{object} has \code{coef} and 
#'  \code{vcov} methods available.
#' 
#' @return A named vector of p-values.
#' @export
pval <- function(object, parm, ...) {
  UseMethod("pval", object)
}

#' @export
pval.default <- function(object, parm) {
  # get parameters
  if (missing(parm)) { 
    parm <- names(coef(object))
  } else if (anyNA(match(parm, names(coef(object))))) {
    stop("Each parameter must match a coefficient name.")
  }

  # Calculate normal approximation p-values
  se <- sqrt(diag(vcov(object))[parm])
  z <- abs(coef(object)[parm]) / se
  return(2 * pnorm(-z))
}

## Define indicator for external data rows
#' Identify rows corresponding to external contact intervals.
#'
#' Specifies indicator variable for rows containing an external contact
#' interval (from the community or environment to a susceptible individual)
#' in model formula.
#'
#' @param ex A binary variable such that \code{as.numeric(ex)} will produce 
#'  a one for each external row and a zero for each internal row
#' @details
#'  If \code{ex} is an indicator variable for external infection, then
#'  \code{exog(ex)} can be added to the model formula to allow separate 
#'  handling of internal and external contact intervals.
#' 
#' @return A factor with two levels
#' @export
ext <- function(x) {
  nx <- as.numeric(x)
  if (!(all(nx %in% c(0, 1)))) stop("External indicator must be binary.")
  return(nx)
}

