#' Methods to work with a `model`.
#' @name methods
#' @description
#' For a given model object query the `initial`,
#' `mu`, log `prior`, `graph` or precision `prec`
#' can be evaluated/retrieved.
#' @param model object to represent a model
#' @return the result of the desired query
#' of the 'cgeneric' model.
#' 'graph' and 'prec' can be either a vector
#' (if optimize = TRUE) or a sparse matrix.
NULL
#> NULL

#' @describeIn methods
#' Retrieve the initial model parameter(s)
#' @export
initial <- function(model) {
  UseMethod("initial")
}
#' @describeIn methods
#' Evaluate the model's mean
#' @param theta numeric vector.
#' For `prior` it can be a numeric matrix,
#' with number of lines equal the size of `theta`
#' and each column as a different case.
#' @export
mu <- function(model, theta) {
  UseMethod("mu")
}
#' @describeIn methods
#' Evaluate the log-prior for a given `theta`
#' @seealso [prior.cgeneric()]
#' @export
prior <- function(model, theta) {
  UseMethod("prior")
}
#' @describeIn methods
#' Retrieve the models' graph
#' @param optimize logical indicating if it is to be
#' returned only the elements and not as a sparse matrix.
#' @export
graph <- function(model, optimize) {
  UseMethod("graph")
}
#' @describeIn methods
#' Retrieve the precision for a given `theta`
#' @export
prec <- function(model, theta, optimize) {
  UseMethod("prec")
}
#' @describeIn methods
#' The default precision method
#' computes the inverse of the variance
#' @param ... additional arguments passed on
#' @export
prec.default <- function(model, ...) {
  v <- vcov(model, ...)
  return(
    forwardsolve(
      backsolve(
        chol(v)
      )
    )
  )
}
#' @describeIn methods
#' The `vcov` method for sparse matrices
#' @param object Matrix supposed to be a
#' sparse precision matrix
setMethod(
  "vcov",
  "Matrix",
  function(object, ...) {
    object <- Matrix::Cholesky(object)
    return(solve(object))
  }
)
#' Define the is.zero method
#' @param x an R object
#' @param ... additional arguments
#' @return logical
#' @export
is.zero <- function(x, ...) {
  UseMethod("is.zero")
}
#' @describeIn is.zero
#' The is.zero.default definition
#' @export
is.zero.default <- function(x, ...) {
  a <- abs(as.numeric(c(x)))
  if(diff(range(a))<(.Machine$double.eps^0.9)) {
    tol <- (.Machine$double.eps^0.9)
  } else {
    tol <- .Machine$double.eps *
      max(sqrt(length(a))) * max(a)
  }
  return(a < tol)
}
#' @describeIn is.zero
#' The is.zero.matrix definition
#' @export
is.zero.matrix <- function(x, ...) {
  stopifnot(inherits(x, "matrix"))
  a <- abs(x)
  if(diff(range(a))<(.Machine$double.eps^0.9)) {
    tol <- (.Machine$double.eps^0.9)
  } else {
    tol <- .Machine$double.eps *
      max(sqrt(length(a))) * max(a)
  }
  return(a < tol)
}

