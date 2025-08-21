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
#' The graph method for 'rgeneric'
#' @param model a `rgeneric` model object
#' @param ... additional arguments
#' @export
graph.rgeneric <- function(model, ...) {
  return(do.call(
    what = "inla.rgeneric.q",
    args = list(rmodel = model,
                cmd = "graph")
    ))
}
#' @describeIn rgeneric-class
#' The precision method for an `rgeneric` object.
#' @param ... additional parameter such as 'theta'
#' If 'theta' is not supplied, initial will be taken.
#' @export
prec.rgeneric <- function(model, ...) {
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
  } else {
    warning("Using the 'default' initial parameter:")
    theta <- initial(model)
    cat(theta, '\n')
  }
  return(do.call(
    what = "inla.rgeneric.q",
    args = list(rmodel = model,
                cmd = "Q",
                theta = theta)
  ))
}
#' @describeIn rgeneric-class
#' The initial method for 'rgeneric'
#' @export
initial.rgeneric <- function(model) {
  return(do.call(
    what = "inla.rgeneric.q",
    args = list(rmodel = model,
                cmd = "initial")
  ))
}
#' @describeIn rgeneric-class
#' The mu method for 'rgeneric'
#' @export
mu.rgeneric <- function(model, theta) {
  if(missing(theta)) {
    theta <- initial(model)
    warning(paste(
      "Using the default initial theta as:",
      format(theta)))
  }
  return(do.call(
    what = "inla.rgeneric.q",
    args = list(rmodel = model,
                cmd = "mu",
                theta = theta)
  ))
}
#' @describeIn rgeneric-class
#' The prior metho for 'rgeneric'
#' @param theta the parameter.
#' @export
prior.rgeneric <- function(model, theta) {
  return(do.call(
    what = "inla.rgeneric.q",
    args = list(rmodel = model,
                cmd = "log.prior",
                theta = theta)
  ))
}
