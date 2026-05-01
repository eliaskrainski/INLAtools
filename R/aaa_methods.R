#' Methods to work with a `model`.
#' @name INLAtools-methods
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

#' @describeIn INLAtools-methods
#' The `vcov` method for sparse matrices
#' @param object Matrix supposed to be a
#' sparse precision matrix
vcov.Matrix <- function(object, ...) {
  object <- Matrix::Cholesky(object)
  return(solve(object))
}
#' Define the is.zero method
#' @param x an R object
#' @param tol numeric to be used as (absolute) tolerance.
#' if missing (default) it will consider `x==0`.
#' @return logical
#' @export
is.zero <- function(x, tol) {
  UseMethod("is.zero")
}
#' @describeIn is.zero
#' The is.zero.default definition
#' @export
is.zero.default <- function(x, tol) {
  if(missing(tol)) {
    return(x==0)
  }
  a <- abs(as.numeric(c(x)))
  return(a < tol)
}
#' @describeIn is.zero
#' The is.zero.matrix definition
#' @export
is.zero.matrix <- function(x, tol) {
  stopifnot(inherits(x, "matrix"))
  return(matrix(is.zero.default(x, tol),
                nrow(x), ncol(x)))
}
#' @describeIn is.zero
#' The is.zero.Matrix definition
#' @export
is.zero.Matrix <- function(x, tol) {
  stopifnot(inherits(x, "Matrix"))
  is.zero(as.matrix(x), tol)
}
