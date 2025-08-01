#' Build an `cgeneric` object for a `generic0` model.
#' See details.
#' @description
#' Build data needed to implement a model whose
#' precision has a conditional precision parameter.
#' This uses the C interface in the 'INLA' package,
#' that can be used as a linear predictor
#' model component with an 'f' term.
#' @param R the structure matrix for the model definition.
#' @param param length two vector with the parameters
#' `a` and `p` for the PC-prior distribution defined from
#'   \deqn{P(\sigma > a) = p}
#' where \eqn{\sigma} can be interpreted as marginal standard
#' deviation of the process if scale = TRUE. See details.
#' @param constr logical indicating if it is to add a
#' sum-to-zero constraint. Default is TRUE.
#' @param scale logical indicating if it is to scale
#' the model. See detais.
#' @param ... arguments (debug,useINLAprecomp,shlib)
#' passed on to [cgeneric()].
#' @details
#' The precision matrix is defined as
#'  \deqn{Q = \tau R}
#' where the structure matrix R is supplied by the user
#' and \eqn{\tau} is the precision parameter.
#' Following Sørbie and Rue (2014), if scale = TRUE
#' the model is scaled so that
#'  \deqn{Q = \tau s R}
#'  where \eqn{s} is the geometric mean of the diagonal
#'  elements of the generalized inverse of \eqn{R}.
#' \deqn{s = \exp{\sum_i \log((R^{-})_{ii})/n}}
#' If the model is scaled, the geometric mean of the
#' marginal variances, the diagonal of \eqn{Q^{-1}}, is one.
#' Therefore, when the model is scaled,
#' \eqn{\tau} is the marginal precision,
#' otherwise \eqn{\tau} is the conditional precision.
#' @references
#' Sigrunn Holbek Sørbye and Håvard Rue (2014).
#' Scaling intrinsic Gaussian Markov random field priors in
#' spatial modelling. Spatial Statistics, vol. 8, p. 39-51.
#' @return a `cgeneric` object, see [cgeneric()].
#' @seealso [prior.cgeneric()]
#' @importFrom methods as
#' @export
#' @examples
#' ## structured precision matrix model definition
#' R <- Matrix(toeplitz(c(2,-1,0,0,0)))
#' R
#' mR <- cgeneric("generic0", R = R,
#'   param = c(1, 0.05), scale = FALSE)
#' graph(mR)
#' prec(mR, theta = 0)
cgeneric_generic0 <-
  function(R,
           param,
           constr = TRUE,
           scale = TRUE,
           ...) {

    stopifnot(param[1]>0)
    if(is.na(param[2])) {
      param[2] = 0.0
    }
    stopifnot(param[2]>=0)
    stopifnot(param[2]<=1)

    dotArgs <- list(...)
    if(is.null(dotArgs$debug)) {
      debug <- FALSE
    } else {
      debug <- dotArgs$debug
    }

    R <- upperPadding(R)
    if(debug) {
      print(str(R))
    }

    n <- as.integer(nrow(R))
    stopifnot(n>0)

    if(scale) {
      Rs <- try(
        INLA::inla.scale.model(
          Q = R,
          constr = list(A = matrix(1, 1, n), e = 0)
      ))
      if(inherits(R, "try-error")) {
        warning("Please install 'INLA'!")
      } else {
        R <- Rs
      }
    }

    if(is.null(dotArgs$useINLAprecomp)) {
      useINLAprecomp <- FALSE
    } else {
      useINLAprecomp <- dotArgs$useINLAprecomp
    }
    if(useINLAprecomp) {
     shlib <- cgeneric_shlib(
       package = "graphpcor",
       useINLAprecomp = TRUE,
       debug = debug)
    } else {
      shlib <- cgeneric_shlib(
        package = "INLAtools",
        useINLAprecomp = FALSE,
        debug = debug)
    }

    the_model <- do.call(
      what = "cgenericBuilder",
      args = list(
        model = "inla_cgeneric_generic0",
        n=as.integer(n),
        param=param,
        Rgraph = R,
        debug = debug,
        shlib = shlib
      )
    )

    if(constr) {
      the_model$f$extraconstr <- list(
        A = matrix(1, 1, n),
        e = 0
      )
    }

    return(the_model)

  }
#' @describeIn cgeneric_generic0
#' The [cgeneric_iid] uses the [cgeneric_generic0]
#' with the structure matrix as the identity.
#' @param n integer required to specify the model size
#' @importFrom Matrix Diagonal
#' @export
cgeneric_iid <-
  function(n,
           param,
           constr = FALSE,
           ...) {
    do.call(
      what = "cgeneric_generic0",
      args = list(
        R = Diagonal(n = n, x = rep(1, n)),
        param = param,
        constr = constr,
        scale = FALSE,
        ...)
    )
  }

