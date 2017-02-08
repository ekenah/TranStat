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
pval <- function(object, parm, ...) {
  UseMethod("pval", object)
}

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

#pairwise <- function(data, cluster) {
#  # generate pairwise data from line data with a cluster variable
#
#}
#
#sus <- function(id) {
#  return(as.factor(id))
#}
