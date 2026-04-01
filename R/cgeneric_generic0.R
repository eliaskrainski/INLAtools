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
#' @importFrom Matrix as.matrix
#' @export
cgeneric_generic0 <-
  function(R,
           param,
           constr = TRUE,
           scale = TRUE,
           ...) {

    stopifnot(param[1]>0)
    if(is.na(param[2])) {
      param[2] <- 0.0
    }
    stopifnot(param[2]>=0)
    stopifnot(param[2]<=1)

    dotArgs <- list(...)
    if(is.null(dotArgs$debug)) {
      dotArgs$debug <- FALSE
    }

    R <- upperPadding(R)
    if(dotArgs$debug) {
      print(str(R))
    }

    n <- as.integer(nrow(R))
    stopifnot(n>0)

    if(scale | constr) {
      ## This work with dense matrices
      ## For big matrices, consider INLA:::inla.scale.model()
      Rd <- as.matrix(R)
      R.svd <- svd(Rd + t(Rd) - diag(diag(Rd)))
      R.svd$ip <- which(R.svd$d>sqrt(.Machine$double.eps))
      if(dotArgs$debug) {
        print(str(list(R.svd=R.svd)))
      }
      stopifnot(length(R.svd$ip)>0)
    }

    if(scale) {
      R.gi <- R.svd$v[, R.svd$ip, drop = FALSE] %*% (
          (1/R.svd$d[R.svd$ip])*t(R.svd$u[, R.svd$ip, drop=FALSE]))
      R <- R * exp(mean(log(diag(R.gi))))
    }
    if(is.null(dotArgs$useINLAprecomp)) {
      dotArgs$useINLAprecomp <- TRUE
    }
    INLAvcheck <- packageCheck("INLA", "25-10-28")
    if(is.na(INLAvcheck) & dotArgs$useINLAprecomp) {
      dotArgs$useINLAprecomp <- FALSE
      warning("INLA version is old. Setting 'useINLAprecomp = FALSE'!")
    }
    dotArgs$shlib <- cgeneric_shlib(
      package = "INLAtools",
      useINLAprecomp = dotArgs$useINLAprecomp,
      debug = dotArgs$debug)

    the_model <- do.call(
      what = "cgenericBuilder",
      args = c(list(
        model = "inla_cgeneric_generic0",
        n=as.integer(n),
        param=param,
        Rgraph = R),
        dotArgs)
    )
## setup extraconstr
    Ae <- dotArgs$extraconstr
    if(!is.null(Ae)) {
      stopifnot(!(names(Ae) %in% c("A", "e")))
      stopifnot(ncol(Ae$A)==n)
      stopifnot(nrow(Ae$A)==length(Ae$e))
    }
    if(constr) {
      the_model$f$constr <- FALSE
      jj <- setdiff(1:n, R.svd$ip)
      Ae$A <- rbind(Ae$A,
                    t(R.svd$u[, jj, drop = FALSE]))
      Ae$e <- c(Ae$e, rep(0, length(jj)))
      qrc <- qr(t(Ae$A))
      c.ok <- which(abs(qrc$qraux)>0)
      Ae$A <- Ae$A[qrc$pivot[c.ok], , drop = FALSE]
      Ae$e <- Ae$e[qrc$pivot[c.ok]]
    }
    the_model$f$extraconstr <- Ae

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

