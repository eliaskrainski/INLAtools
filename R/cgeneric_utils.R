#' `cgeneric_get` is an internal function used by
#' `graph`, `pred`, `initial`, `mu` or `prior`
#' methods for `cgeneric`.
#' @description
#' The `generic_get` retrieve a model property specified by
#' `cmd` on an `cgeneric` object.
#' The functions listed below are for each `cmd` case.
#' @param model a `cgeneric` object.
#' @param cmd an string to specify which model element to get
#' @param theta numeric vector with the model parameters.
#' If missing, the [initial()] will be used.
#' @param optimize logical. If missing or FALSE,
#' the graph and precision are as a sparse matrix.
#' If TRUE, graph only return the row/col indexes and
#' precision return only the elements as a vector.
#' @useDynLib INLAtools
#' @return depends on `cmd`
#' @seealso check the examples in [cgeneric_generic0()]
cgeneric_get <- function(model,
                         cmd = c("graph", "Q", "initial", "mu", "log_prior"),
                         theta,
                         optimize = TRUE) {

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

  ##  print(c(cmd = cmd))

  cgdata <- model$f$cgeneric$data
  stopifnot(!is.null(cgdata))
  stopifnot(!is.null(cgdata$ints))
  stopifnot(!is.null(cgdata$characters))

  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

    if(missing(theta)) {
      if(cmd %in% c("Q", "log_prior")) {
        stop("Please provide 'theta'!")
      } else {
        theta <- NULL
        ntheta = 0L
      }
    } else {
      storage.mode(theta) <- "double"
      if(inherits(theta, "matrix")) {
        ntheta <- as.integer(ncol(theta))
      } else {
        ntheta <- 1L
      }
    }

    if(length(cmd) == 1) {
      ret <- .Call(
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
      )

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
        ret <- Matrix::sparseMatrix(
          i = ij[[1]] + 1L,
          j = ij[[2]] + 1L,
          x = ret,
          symmetric = TRUE,
          repr = "T"
        )
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

    if(any(cmd == "graph")) {
      ret$graph <-
        Matrix::sparseMatrix(
          i = ret$graph[[1]] + 1L,
          j = ret$graph[[2]] + 1L,
          x = rep(1, length(ret$graph[[1]])),
          symmetric = TRUE,
          repr = "T"
        )
    }

    if(any(cmd == "Q")) {
      if(any(cmd == "graph")) {
        ij <- ret$graph
      } else {
        ij <- .Call(
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
        )
        ij <- Matrix::sparseMatrix(
          i = ij[[1]] + 1L,
          j = ij[[2]] + 1L,
          symmetric = TRUE,
          repr = "T"
        )
      }
      ij@x <- ret$Q
      ret$Q <- ij
    }

  return(ret)

}
#' @describeIn cgeneric_get
#' Retrive the initial parameter(s) of an `cgeneric` model.
#' @export
initial.cgeneric <- function(model) {
  cgeneric_get(model, "initial")
}
#' @describeIn cgeneric_get
#' Evaluate the mean for an `cgeneric` model.
#' @export
mu.cgeneric <- function(model, theta) {
  cgeneric_get(model, "mu", theta = theta)
}
#' @describeIn cgeneric_get
#' Evaluate the prior for an `cgeneric` model
#' @return numeric scalar (if numeric vector is provided
#' for theta) or vector (if numeric matrix is provided
#' for theta).
#' @export
#' @seealso [cgeneric_generic0()]
#' @examples
#' old.par <- par(no.readonly = TRUE)
#'
#' ## Setting the prior parameters
#' prior.par <- c(1, 0.5) # P(sigma > 1) = 0.5
#' cmodel <- cgeneric(
#'   model = "iid", n = 10,
#'   param = c(prior.par), useINLAprecomp = FALSE)
#'
#' ## prior summaries: sigma and log-precision
#' (lamb <- -log(prior.par[2])/prior.par[1])
#' (smedian <- qexp(0.5, lamb))
#' (smean <- 1/lamb)
#'
#' ## mode: at the minimum of - log-prior
#' (lpmode <- optimize(function(x)
#'   -prior(cmodel, theta = x),
#'   c(-10, 30))$minimum)
#' ## mean: integral of x*f(x)dx
#' (lpmean <- integrate(function(x)
#'   exp(prior(cmodel, theta = matrix(x, 1)))*x,
#'   -10, 30)$value)
#'
#' ## prior visualization: log(precision) and sigma
#' par(mfrow = c(1, 2))
#' plot(function(x)
#'  exp(prior(cmodel, theta = matrix(x, nrow=1))),
#'   -3, 3, n = 601, xlab = "log-precision",
#'   ylab = "density")
#' abline(v = lpmode, lwd = 3, col = 2)
#' rug(-2*log(smedian), lwd = 3, col = 3)
#' rug(lpmean, lwd = 3, col = 4)
#' plot(function(x)
#'  exp(prior(cmodel,
#'   theta = matrix(
#'     -2*log(x),
#'     nrow = 1))+log(2)-log(x)),
#'   1/100, 10, n = 1000,
#'   xlab = expression(sigma),
#'   ylab = "density")
#' plot(function(x) dexp(x, lamb),
#'    1/100, 10, n = 1000,
#'    add = TRUE, lty = 2, col = 2)
#' rug(smedian, lwd = 3, col = 3)
#' rug(smean, lwd = 3, col = 4)
#' par(old.par)
prior.cgeneric <- function(model, theta) {
  return(cgeneric_get(model = model,
                      cmd = "log_prior",
                      theta = theta))
}
#' @describeIn cgeneric_get
#' Retrieve the graph of an `cgeneric` object
#' @param optimize logical indicating if it is to be
#' returned only the elements and not as a sparse matrix.
#' @export
graph.cgeneric <- function(model, optimize) {
  if(missing(optimize)) {
    optimize <- FALSE
  }
  stopifnot(is.logical(optimize))
  return(cgeneric_get(
    model, "graph",
    optimize = optimize))
}
#' @describeIn cgeneric_get
#' Retrieve the precision of an `cgeneric` object
#' @export
prec.cgeneric <- function(model, theta, optimize) {
  if(missing(optimize)) {
    optimize <- FALSE
  }
  stopifnot(is.logical(optimize))
  cgeneric_get(model, cmd = "Q", theta = theta, optimize = optimize)
}
