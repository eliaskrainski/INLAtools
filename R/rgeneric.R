#' @rdname rgeneric-class
#' @param model an object used to define the model.
#' See the 'rgeneric' vignette from the INLA package.
#' @param debug logical indicating debug state.
#' @param n integer with the dimension of the model
#' @param ... additional arguments to be used internally
#' for the model, for example, additional data.
#' @returns `rgeneric`/ `inla.rgeneric` object.
#' @export
rgeneric <- function(model,
                     n,
                     debug = FALSE,
                     ...) {
  UseMethod("rgeneric")
}
#' The rgeneric default method.
#' @rdname rgeneric-class
#' @export
rgeneric.default <- function(model,
                             n,
                             debug = FALSE,
                             ...) {
  ## as INLA::inla.rgeneric.define(..., compile = TRUE, optimize = TRUE)
  rmodel <- structure(
    list(
      f = list(
        model = "rgeneric",
        n = n,
        rgeneric = structure(
          list(
            definition =
              compiler::cmpfun(
                model,
                options = list(optimize = 3L)),
            debug = debug,
            optimize = TRUE
          ),
          # inla.rgeneric is needed to support INLA before August 2025
          class = c("inla.rgeneric.f", "inla.rgeneric")
        )
      )
    ),
    class = c("rgeneric", "inla.rgeneric")
  )

  return(rmodel)
}

#' @describeIn rgeneric-class Returns the model object unchanged.
#' @export
rgeneric.rgeneric <- function(model, ...) {
  return(model)
}

#' @describeIn rgeneric-class Converts a regular `inla.rgeneric` object to `rgeneric`.
#' @export
rgeneric.inla.rgeneric <- function(model, ...) {
  class(model) <- c("rgeneric", class(model))
  return(model)
}

#' @describeIn rgeneric-class
#' Print the rgeneric object
#' @param x a rgeneric object
#' @param ... not used
#' @export
print.rgeneric <- function(x, ...) {
  cat("rgeneric: ", x$f$rgeneric$model, ", n = ",
      x$f$rgeneric$n, "\n", sep = "")
}
#' @describeIn rgeneric-class
#' A summary for a rgeneric object
#' @param object a rgeneric object
#' @param ... not used
#' @export
summary.rgeneric <- function(object, ...) {
  g <- graph(object)
  cat("n = ", object$f$rgeneric$n, ", graph with",
      length(g@x), "non-zeros\n", sep = "")
}
#' @describeIn rgeneric-class
#' A plot for a rgeneric object
#' @param y not used
#' @importFrom graphics image
#' @export
plot.rgeneric <- function(x, y, ...) {
  g <- graph(x)
  image(g)
}
