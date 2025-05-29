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
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param libpath string, default is NULL, with the path to the shared object.
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
#' @useDynLib INLAtools
#' @importFrom methods as
#' @examples
#' ## structured precision matrix model definition
#' R <- Matrix(toeplitz(c(2,-1,0,0,0)))
#' mR <- cgeneric("generic0", R = R,
#'   scale = FALSE, param = c(1, 0.05),
#'   useINLAprecomp = FALSE)
#' graph(mR)
#' prec(mR, theta = 0)
cgeneric_generic0 <-
  function(R,
           param,
           constr = TRUE,
           scale = TRUE,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    if(is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("INLAtools")
      } else {
        libpath <- system.file("libs", package = "INLAtools")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "x64/INLAtools.dll")
        } else {
          libpath <- file.path(libpath, "INLAtools.so")
        }
      }
    }

    stopifnot(param[1]>0)
    if(is.na(param[2])) {
      param[2] = 0.0
    }
    stopifnot(param[2]>=0)
    stopifnot(param[2]<=1)

    R <- Sparse(R)

    n <- as.integer(nrow(R))
    stopifnot(n>0)

    idx <- which(R@i <= R@j)

    if(debug) {
      print(str(list(
        ii = R@i,
        jj = R@j,
        idx = idx
      )))
    }

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

    ord <- order(R@i[idx])
    nnz <- length(idx)
    cmodel = "inla_cgeneric_generic0"

    the_model <- list(
      f = list(
        model = "cgeneric",
        n = n,
        cgeneric = list(
          model = cmodel,
          shlib = libpath,
          n = as.integer(n),
          debug = as.logical(debug),
          data = list(
            ints = list(
              n = as.integer(n),
              debug = as.integer(debug)
            ),
            doubles = list(
              param = param
            ),
            characters = list(
              model = cmodel,
              shlib = libpath
            ),
            matrices = list(
            ),
            smatrices = list(
              Rgraph = c(
                n, n, nnz,
                R@i[idx][ord],
                R@j[idx][ord],
                R@x[idx][ord]
              )
            )
          )
        )
      )
    )

    class(the_model) <- "cgeneric"
    class(the_model$f$cgeneric) <- "cgeneric"

    if(constr) {
      the_model$f$extraconstr <- list(
        A = matrix(1, 1, n),
        e = 0
      )
    }

    return(the_model)

  }
#' The [cgeneric_iid] uses the [cgeneric_generic0]
#' with the structure matrix as the identity.
#' @param n integer required to specify the model size
#' @inheritParams cgeneric_generic0
#' @importFrom Matrix Diagonal
cgeneric_iid <-
  function(n,
           param,
           constr = FALSE,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {
    cgeneric_generic0(
      R = Diagonal(n = n,
                   x = rep(1, n)),
      param = param,
      constr = constr,
      scale = FALSE,
      debug = debug,
      useINLAprecomp = useINLAprecomp,
      libpath = libpath
    )
}

