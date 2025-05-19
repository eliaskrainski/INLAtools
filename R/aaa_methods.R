#' Methods to work with a `model`.
#' @name methods
#' @description
#' For a given model object the `intial`,
#' `mu`, log `prior`, `graph` or precision `prec`
#' can be evaluated/retrieved, each one
#' through a method function.
#' @param model object to represent a model
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
