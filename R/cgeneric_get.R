#' `cgeneric_get` is an internal function used to query
#' `graph`, `Q`, `initial`, `mu` or `log_prior` from a
#' `cgeneric` model.
#' @description
#' The `generic_get` retrieve a model property specified by
#' `cmd` on an `cgeneric` object.
#' The functions listed below are for each `cmd` case.
#' @param model a `cgeneric` object.
#' @param cmd an string to specify which model element to get
#' @param theta numeric vector with the model parameters.
#' If missing, the `initial` will be used.
#' @param optimize logical. If missing or FALSE,
#' the graph and precision are as a sparse matrix.
#' If TRUE, graph only return the row/col indexes and
#' precision return only the elements as a vector.
#' @return depends on `cmd`
#' @seealso check the examples in [cgeneric_generic0()]
cgeneric_get <- function(model,
                         cmd = c("graph", "Q", "initial", "mu", "log_prior"),
                         theta,
                         optimize = TRUE) {

  stopifnot(model$f$n>0)
  n <- model$f$n

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

  ##  print(c(cmd = cmd))

  cgdata <- model$f$cgeneric$data
  stopifnot(!is.null(cgdata))
  stopifnot(!is.null(cgdata$ints))
  stopifnot(!is.null(cgdata$characters))

  ## from Version 0.0.3.902 (and src/cgeneric_get.c)
  ## split each sparse matrix into list elements:
  ## nr, nc, m, i, j, x
  nsm <- length(cgdata$smatrices)
  if(nsm>0) {
    for(i in 1:nsm) {
      smi <- cgdata$smatrices[[i]]
      mi <- as.integer(smi[3])
      cgdata$smatrices[[i]] <- list(
        nr = as.integer(smi[1]),
        nc = as.integer(smi[2]),
        m = mi,
        i = as.integer(smi[3+1:mi]),
        j = as.integer(smi[3+mi+1:mi]),
        x = smi[3+2*mi+1:mi]
      )
    }
  }

  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

  initheta <- try(.Call(
    "inla_cgeneric_element_get",
    "initial",
    NULL,
    as.integer(1),
    cgdata$ints,
    cgdata$doubles,
    cgdata$characters,
    cgdata$matrices,
    cgdata$smatrices,
    PACKAGE = "INLAtools"
  ), silent = TRUE)
  if(inherits(initheta, "try-error")) {
    cat("Problem with the current `shlib`:")
    print(attr(initheta, "condition")$message)
    ## workaround: swap   shlib
    ## Note: cgeneric Kronecker not considered
    sshlib <- strsplit(cgdata$characters$shlib, "/")[[1]]
    lpkg <- utils::tail(sshlib, 2)[1]
    ish <- cgeneric_shlib_path(
      package = lpkg,
      useINLAprecomp = FALSE)
    psh <- cgeneric_shlib_path(
      package = lpkg,
      useINLAprecomp = FALSE)
    if(cgdata$characters$shlib == ish) {
      cgdata$characters$shlib <- psh
    } else {
      cgdata$characters$shlib <- ish
    }
    warning(paste("Changed `shlib` to\n",
                  cgdata$characters$shlib))
    initheta <- try(.Call(
      "inla_cgeneric_element_get",
      "initial",
      NULL,
      as.integer(1),
      cgdata$ints,
      cgdata$doubles,
      cgdata$characters,
      cgdata$matrices,
      cgdata$smatrices,
      PACKAGE = "INLAtools"
    ), silent = TRUE)
    if(inherits(initheta, "try-error")) {
      print(attr(initheta, "condition")$message)
      stop('Error trying to get "initial"!')
    }
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
        warning('NULL "theta", using "initial"!')
        theta <- initheta
      }
    }
    theta <- as.double(theta)
    stopifnot((length(theta)%%theta.size)==0)
    ntheta <- floor(length(theta)/length(initheta))
  } else {
    theta <- NULL
    ntheta <- 0L
  }

  if(any(cmd %in% c("graph", "Q")) & (!optimize)) {
    ij2Q <- function(ij, x = NULL) {
      if(is.null(x)) {
        x <- rep(1L, length(ij[[1]]))
      }
      Sparse(Matrix::sparseMatrix(
        i = ij[[1]] + 1L,
        j = ij[[2]] + 1L,
        x = x,
        symmetric = TRUE,
        repr = "T",
        dims = c(n, n)
      ))
    }
  }

  if(length(cmd) == 1) {
    ret <- try(.Call(
      "inla_cgeneric_element_get",
      cmd,
      theta,
      as.integer(ntheta),
      cgdata$ints,
      cgdata$doubles,
      cgdata$characters,
      cgdata$matrices,
      cgdata$smatrices,
      PACKAGE = "INLAtools"
    ), silent = TRUE)
    if(inherits(ret, "try-error")) {
      print(attr(ret, "condition")$message)
      stop('Error trying to get "', cmd, '"!')
    }

    if((cmd %in% c("graph", "Q")) && (!optimize)) {
      if(cmd == "graph") {
        ij <- ret
        ret <- rep(1, length(ij[[1]]))
      } else {
        ij <- .Call(
          "inla_cgeneric_element_get",
          "graph",
          NULL,
          as.integer(ntheta),
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "INLAtools"
        )
      }
      ret <- ij2Q(ij, ret)
    }
    return(ret)
  }

  names(cmd) <- cmd
  ret <-
    lapply(
      cmd, function(x) {
        .Call(
          "inla_cgeneric_element_get",
          x,
          theta,
          as.integer(ntheta),
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "INLAtools"
        )
      }
    )
  if(optimize) {
    return(ret)
  }

  if(any(cmd == "Q")) {
    if(any(cmd == "graph")) {
      ret$Q <- ret$graph <- ij2Q(ret$graph, ret$Q)
      ret$graph@x <- rep(1L, length(ret$graph@x))
    } else {
      ret$Q <- ij2Q(
        .Call(
          "inla_cgeneric_element_get",
          "graph",
          theta,
          as.integer(ntheta),
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "INLAtools"
        ),
        x = ret$Q
      )
    }
  }

  return(ret)

}
#' @describeIn cgeneric_get
#' Retrive the initial parameter(s) of an `cgeneric` model.
#' @export
cgeneric_initial <- function(model) {
  cgeneric_get(model, "initial")
}
#' @describeIn cgeneric_get
#' Evaluate the mean for an `cgeneric` model.
#' @export
cgeneric_mu <- function(model, theta) {
  cgeneric_get(model, "mu", theta = theta)
}
#' @describeIn cgeneric_get
#' Retrieve the graph of an `cgeneric` object
#' @param optimize logical indicating if it is to be
#' returned only the elements and not as a sparse matrix.
#' @export
cgeneric_graph <- function(model, optimize) {
  if(missing(optimize)) {
    optimize <- FALSE
  }
  return(cgeneric_get(
    model, "graph",
    optimize = optimize))
}
#' @describeIn cgeneric_get
#' Retrieve the precision of an `cgeneric` object
#' @export
cgeneric_Q <- function(model, theta, optimize) {
  if(missing(optimize)) {
    optimize <- FALSE
  }
  cgeneric_get(model,
               cmd = "Q",
               theta = theta,
               optimize = optimize)
}
#' @describeIn cgeneric_get
#' Evaluate the prior for an `cgeneric` model
#' @return numeric scalar (if numeric vector is provided
#' for theta) or vector (if numeric matrix is provided
#' for theta).
#' @export
#' @example demo/prior.R
cgeneric_prior <- function(model, theta) {
  return(cgeneric_get(model = model,
                      cmd = "log_prior",
                      theta = theta))
}
