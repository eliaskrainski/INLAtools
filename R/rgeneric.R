#' @rdname rgeneric-class
#' @param model an object used to define the model.
#' Its class will define which method is considered.
#' @param debug logical indicating debug state.
#' @param compile logical indicating to compile the model.
#' @param optimize logical indicating if only the elements
#' of the precision matrix are returned.
#' @param ... additional arguments to be used internally
#' for the model, for example, additional data.
#' @returns `rgeneric`/ `inla.rgeneric` object.
#' @export
rgeneric <- function(model,
                     debug = FALSE,
                     compile = TRUE,
                     optimize = TRUE,
                     ...) {
  UseMethod("rgeneric")
}
#' The rgeneric default method.
#' @rdname rgeneric-class
#' @param model For the default method, the model defined as a function.
#' See the 'rgeneric' vignette from the INLA package.
#' @export
rgeneric.default <- function(model,
                             debug = FALSE,
                             compile = TRUE,
                             optimize = TRUE,
                             ...) {
  ## it uses INLA::inla.rgeneric.define()

  rmodel <- INLA::inla.rgeneric.define(
    model = model,
    debug = debug,
    compile = compile,
    optimize = optimize,
    ...
  )

  class(rmodel) <- c("rgeneric", class(rmodel))
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
