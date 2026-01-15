#' `rgeneric_get` is an internal function used by
#' `graph`, `pred`, `initial`, `mu` or `prior`
#' methods for `rgeneric`.
#' @description
#' The `generic_get` retrieve a model property specified by
#' `cmd` on an `rgeneric` object.
#' The functions listed below are for each `cmd` case.
#' @param model a `rgeneric` object.
#' @param cmd an string to specify which model element to get
#' @param theta numeric vector with the model parameters.
#' If missing, the [initial()] will be used.
#' @param optimize logical. If missing or FALSE,
#' the graph and precision are as a sparse matrix.
#' If TRUE, graph only return the row/col indexes and
#' precision return only the elements as a vector.
#' @return depends on `cmd`
rgeneric_get <- function(model,
                         cmd = c("graph", "Q", "initial", "mu", "log_prior"),
                         theta,
                         optimize = TRUE) {

  ## internal function to
  ## re-implement INLA::inla.rgeneric.q
  ## to be accessed by rgeneric INLAtools methods

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)
  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

  func <- model$f$rgeneric$definition
  initheta <- try(do.call(
    what = func,
    args = list(cmd = "initial", theta = NULL)
  ), silent = TRUE)

  if(inherits(initheta, "try-error")) {
    stop('Error trying to get "initial"!')
  }
  if((length(cmd)==1) && (cmd=="initial")) {
    return(initheta)
  }
  theta.size <- length(initheta)

  needtheta <- any(cmd %in% c("Q", "log_prior"))
  if(needtheta & (theta.size>0)) {
    if(missing(theta)) {
      warning('missing "theta", using "initial"!')
      theta <- initheta
    } else {
      if(is.null(theta)) {
        warning('missing "theta", using "initial"!')
        theta <- initheta
      }
    }
    theta <- as.double(theta)
    ntheta <- floor(length(theta)/length(initheta))

    ## if more than one theta is given
    if((ntheta>1) && (length(cmd)==1) && (cmd=="log_prior")) {
      theta <- matrix(theta, nrow = length(initheta))
      return(sapply(1:ncol(theta), function(j) {
        rgeneric_get(model, cmd = "initial", theta = theta[, j])
      }))
    }

  } else {
    theta <- NULL
    ntheta <- 0L
  }

  if(length(cmd) == 1) {
    ret <- try(do.call(
      what = func,
      args = list(cmd = cmd, theta = theta)
    ), silent = TRUE)
    if(inherits(ret, "try-error")) {
      stop('Error trying to get "', cmd, '"!')
    }
    if(cmd %in% c("graph", "Q"))
      ret <- Sparse(ret)
    if(optimize) {
      if(cmd == "graph")
        ret <- list(i = as.integer(ret@i),
                    j = as.integer(ret@j))
      if(cmd == "Q")
        ret <- ret@x
    }
  } else {
    names(cmd) <- cmd
    ret <- lapply(cmd, function(x) {
      try(do.call(what = func,
                  args = list(cmd = cmd, theta = theta)),
          silent = TRUE)
    })
    if(optimize) {
      if(any(cmd == "graph"))
        if(!inherits(ret$graph, "try-error"))
          ret$graph <- list(
            i = ret$graph@i,
            j = ret$graph@j)
      if(any(cmd == "Q"))
        if(!inherits(ret$Q, "try-error"))
          ret$Q <- ret$Q@x
    }
  }

  return(ret)

}
#' @describeIn rgeneric_get
#' Retrive the initial parameter(s) of an `rgeneric` model.
#' @export
initial.rgeneric <- function(model) {
  rgeneric_get(model, "initial")
}
#' @describeIn rgeneric_get
#' Evaluate the mean for an `rgeneric` model.
#' @export
mu.rgeneric <- function(model, theta) {
  rgeneric_get(model, "mu", theta = theta)
}
#' @describeIn rgeneric_get
#' Retrieve the graph of an `rgeneric` object
#' @param optimize logical indicating if it is to be
#' returned only the elements and not as a sparse matrix.
#' @export
graph.rgeneric <- function(model, optimize) {
  if(missing(optimize)) {
    optimize <- FALSE
  }
  return(rgeneric_get(
    model, "graph",
    optimize = optimize))
}
#' @describeIn rgeneric_get
#' Retrieve the precision of an `rgeneric` object
#' @export
prec.rgeneric <- function(model, theta, optimize) {
  if(missing(theta)) {
    warning('missing "theta", using "initial"!')
    theta <- initial(model)
  }
  if(missing(optimize)) {
    optimize <- FALSE
  }
  rgeneric_get(model,
               cmd = "Q",
               theta = theta,
               optimize = optimize)
}
#' @describeIn rgeneric_get
#' Evaluate the prior for an `rgeneric` model
#' @return numeric scalar (if numeric vector is provided
#' for theta) or vector (if numeric matrix is provided
#' for theta).
#' @export
#' @example demo/prior.R
prior.rgeneric <- function(model, theta) {
  return(rgeneric_get(model = model,
                      cmd = "log_prior",
                      theta = theta))
}
