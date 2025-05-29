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
