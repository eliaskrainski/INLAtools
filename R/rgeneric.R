#' @rdname rgeneric-class
#' @inherit cgeneric-class description
#' @param model an object used to define the model.
#' See the 'rgeneric' vignette from the INLA package.
#' @param debug logical indicating debug state.
#' @param n integer with the dimension of the model
#' @param ... additional arguments to be used internally
#' for the model, for example, additional data.
#' @inherit cgeneric details
#' @note
#' Recommended for prototyping, whereas
#' `cgeneric` is recommended for production.
#' @returns `rgeneric`/ `inla.rgeneric` object.
#' @export
rgeneric <- function(model,
                     n,
                     debug = FALSE,
                     ...) {
  UseMethod("rgeneric")
}
#' @describeIn rgeneric
#' The rgeneric default method.
#' @export
rgeneric.default <- function(model,
                             n,
                             debug = FALSE,
                             ...) {
    dArgs <- list(...)
    stopifnot(!any(names(dArgs)==""))

    env <- if (length(dArgs) > 0) as.environment(dArgs) else new.env()
    parent.env(env) <- .GlobalEnv
    environment(model) <- env

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

#' @describeIn rgeneric Returns the model object unchanged.
#' @export
rgeneric.rgeneric <- function(model, ...) {
  return(model)
}

#' @describeIn rgeneric Check and converts a regular `inla.rgeneric` object to `rgeneric`.
#' @export
rgeneric.inla.rgeneric <- function(model, ...) {
  stopifnot(!c("f") %in% names(model))
  stopifnot(!all(c("model", "n", "rgeneric") %in%
                   names(model$f)))
  stopifnot(!all(c("definition", "debug", "optimize") %in%
                   names(model$f$rgeneric)))
  class(model) <- c("rgeneric", class(model))
  return(model)
}

#' @describeIn rgeneric Print the rgeneric object
#' @param x a rgeneric object
#' @param ... not used
#' @export
print.rgeneric <- function(x, ...) {
  cat("rgeneric: ", x$f$rgeneric$model, ", n = ",
      x$f$rgeneric$n, "\n", sep = "")
}
#' @describeIn rgeneric A summary for a rgeneric object
#' @param object a rgeneric object
#' @param ... not used
#' @export
summary.rgeneric <- function(object, ...) {
  g <- rgeneric_get(object, "graph")
  cat("n = ", object$f$rgeneric$n, ", graph with",
      length(g@x), "non-zeros\n", sep = "")
}
#' @describeIn rgeneric A plot for a rgeneric object
#' @param y not used
#' @export
plot.rgeneric <- function(x, y, ...) {
  g <- rgeneric_get(x, "graph")
  image(g)
}
