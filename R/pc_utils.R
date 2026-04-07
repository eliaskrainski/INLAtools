#' Internal functions to check PC-prior parameters.
#' @name pc-utils
NULL
#> NULL

#' @describeIn pc-utils
#' Check the PC-prior arguments.
#' @param npars number of parameters.
#' @param reference numeric vector to set the reference
#' for each parameter for its PC-prior.
#' @param probability numeric vector with to
#' set the probability statement of the PC prior for each
#' parameter. For `sigma` probability statement is
#' P(sigma > `reference`) = p whereas for `range` it is
#' P(range < `reference`). If  NA, 0 or 1, the corresponding
#' `reference` will be used as fixed. If missing, all the
#' parameters considered as known (fixed) and equal the
#' corresponding reference value.
#' @export
pcParamCheck <- function(npars,
                         reference,
                         probability) {
  if(missing(reference)) {
    reference <- rep(1, npars)
  }
  if(length(reference)<1) {
    reference <-
      rep(reference, npars)[1:npars]
  }
  if(missing(probability)) {
    probability <- rep(0, npars)
  }
  if(length(probability)<1) {
    probability <-
      rep(probability, npars)[1:npars]
  }
  stopifnot(length(reference) == npars)
  stopifnot(length(probability) == npars)
  stopifnot(all(reference>0))
  probability[is.na(probability)] <- 0
  stopifnot(all(probability>=0.0))
  stopifnot(all(probability<=1.0))
  fixed <- is.zero(probability) |
    is.zero(1-probability)
  return(list(
    reference = reference,
    probability = probability,
    fixed = fixed
  ))
}
#' @describeIn pc-utils
#' Penalized Complexity (PC) prior for the log of the practical range.
#' @param lrange numeric with the log of the (practical) range
#' @param lam numeric with the prior parameter
#' @param d integer to specify the domain dimention
#' @param log logical indicating if the density
#' is to be returned in the log scale
#' @export
pclrange <- function(lrange, lam, d = 2, log = FALSE) {
  dh <- 0.5 * d
  out <- log(lam * dh) -dh * lrange - lam * exp(-dh * lrange)
  if(log)
    return(out)
  return(exp(out))
}
#' @describeIn pc-utils
#' Penalized Complexity (PC) prior for the practical range.
#' @param range numeric with the of the (practical) range.
#' @export
#' @examples
#' # P(range < 2.0) = 0.1
#'  lam <- -log(0.1) * 2.0
#'  plot(function(x) pcrange(x, lam), 1/100, 10, n = 100)
pcrange <- function(range, lam, d = 2, log = FALSE) {
  pclrange(log(range), lam, d = d, log = log)/range
}
