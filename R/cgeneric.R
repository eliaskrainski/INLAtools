#' @title Defines a GMRF model to be used with the
#' C interface for `INLA` as a latent model.
#' @description
#' This prepare data for the C type to organize data needed
#' for building latent models which are characterized
#' from given model parameters \eqn{\theta} and the
#' the following model elements.
#'  *  `graph` to define the non-zero precision matrix pattern.
#'  only the upper triangle including the diagonal is needed.
#'  The order should be by line.
#'  * `Q` vector where the
#'     * first element (N) is the size of the matrix,
#'     * second element (M) is the number of non-zero
#'     elements in the upper part (including) diagonal
#'     * the remaining (M) elements are the actual
#'     precision (upper triangle plus diagonal) elements
#'     whose order shall follow the graph definition.
#'  * `mu` the mean vector,
#'  * `initial` vector with
#'    * first element as the number of the parameters in the model
#'    * remaining elements should be the initials for the model parameters.
#'  * `log.norm.const` log of the normalizing constant.
#'  * `log.prior` log of the prior for the model parameters.
#'
#' See details in [INLA::cgeneric()]
#' @param model object class for what a `cgeneric` method exists.
#' if it is a character, a specific function will be called,
#' for example cgeneric("iid", ...") calls cgeneric_iid(...),
#' see [cgeneric_iid()] and [cgeneric_generic0()].
#' @param ... additional arguments passed on to methods
#' @return  named list of `cgeneric` class containing
#'  the named list `f` that contain `model` (a character
#'  always equal to `cgeneric`), `n` (integer)
#'  and `cgeneric` as a named list that contains the
#'  data needed to define the model. Each element on
#'  ...$f$cgeneric is also a named list containing
#'  `ints`, `doubles`, `characters`, `matrices`
#'  and `smatrices`.
#' @seealso [INLA::cgeneric()] and [methods()]
#' @export
cgeneric <- function(model, ...) {
  UseMethod("cgeneric")
}
#' @describeIn cgeneric
#' This calls [INLA::inla.cgeneric.define()]
#' @param model object class for what a `cgeneric` method exists.
#' E.g., if it is a character, a specific function will be called:
#'  cgeneric("iid", ...") calls cgeneric_iid(...)
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param pkg character to give the package name, used to build the
#' libpath (if this is not provided).
#' @param libpath string, default is NULL, with the path to the shared object.
#' @param ... additional arguments passed on to methods
#' @importFrom methods is
#' @export
cgeneric.default <- function(model,
                             debug = FALSE,
                             useINLAprecomp = TRUE,
                             pkg = NULL,
                             libpath = NULL,
                             ...) {
  ## do what INLA::inla.cgeneric.define() does
  if(is.null(libpath)) {
    if (useINLAprecomp) {
	    stopifnot(!is.null(pkg))
      shlib <- INLA::inla.external.lib(pkg)
    } else {
      libpath <- system.file("libs", package = pkg)
      if (Sys.info()["sysname"] == "Windows") {
        shlib <- file.path(libpath,
			   paste0("x64/", pkg, ".dll"))
      } else {
        shlib <- file.path(libpath,
			   paste0(pkg, ".so"))
      }
    }
  } else {
    shlib <- libpath
  }

  args <- c(model = as.character(model),
	    debug = as.integer(debug),
	    shlib = as.character(shlib),
	    list(...))
  nargs <- names(args)
  if(any(nargs) == "")
    stop("Please name the arguments!")
  if(!any(nargs == "n"))
     stop("Please provid 'n'!")
  storage.mode(args$n) <- "integer"
  stopifnot(args$n>=1)

  nM <- length(iM <- which(sapply(args, function(x) is(x, "Matrix"))))
  nm <- length(im <- which(sapply(args, is.matrix)))
  ni <- length(ii <- setdiff(which(sapply(args, is.integer)), im))
  nd <- length(id <- setdiff(which(sapply(args, is.double)), im))
  nc <- length(ic <- which(sapply(args, is.character)))
  stopifnot(ni>1)
  stopifnot(nc>1)

  cmodel <- list(f = list(model = "cgeneric", n=args$n))
  cmodel$cgeneric <- args[c("model", "shlib", "n", "debug")]
  cmodel$data <- vector("list", 5L)
  names(cmodel$data) <- c("ints", "doubles", "characters",
			  "matrices", "smatrices")
  cmodel$data$ints <- args[ii]
  if(nd>0) cmodel$data$doubles <- args[id]
  cmodel$data$characters <- args[ic]
  if(nm>0) cmodel$data$matrices <- args[im]
  if(nM>0) {
	  cmodel$data$smatrices <- args[iM]
	  for(i in 1:nM)
		  cmodel$data$smatrices[[i]] <-
			  Sparse(cmodel$data$smatrices[[i]])
  }
  class(cmodel) <- c("cgeneric",
		     "inla.cgeneric") ## this is needed in INLA::f()
  class(cmodel$f$cgeneric) <- class(cmodel)
  return(cmodel)
}
#' @describeIn cgeneric
#' Method for when `model` is a character.
#' E.g. cgeneric(model = "generic0")
#' calls [cgeneric_generic0]
#' @importFrom methods existsFunction
#' @export
cgeneric.character <- function(model, ...) {
  fn <- paste0("cgeneric_", model)
  if(!existsFunction(fn)) {
    fn <- "cgeneric.default"
  }
  return(do.call(what = fn,
                 args = list(...)))
}
#' Draw samples from hyperparameters of a `cgeneric`
#' model component from an `inla` output, like
#' `inla::inla.iidkd.sample()`.
#' @param n integer as the sample size.
#' @param result an `inla` output.
#' @param name character with the name of the model
#' component in the set of random effects.
#' @param model a `cgeneric` model
#' @param from.theta a function to convert from
#' theta to the desired output for each sample.
#' @param simplify logical (see ?sapply).
#' @return matrix (if n>1 and length(from.theta)>1)
#' or numeric vector otherwise.
#' @seealso [prior.cgeneric()]
#' @export
inla.cgeneric.sample <- function(n = 1e4, result, name,
                                 model, from.theta,
                                 simplify = FALSE) {
  stopifnot(!missing(result))
  stopifnot(inherits(result, "inla"))
  stopifnot(!missing(name))
  stopifnot(!missing(model))
  stopifnot(n > 0)
  idx <- grep(name, names(result$mode$theta))
  xx <- INLA::inla.hyperpar.sample(
    n = n,
    result = result,
    intern = TRUE)[, idx, drop = FALSE]
  if(inherits(model, "cgeneric")) {
    result <- sapply(1:n, function(i)
      cgeneric_get(model = model,
                   cmd = "Q",
                   theta = xx[i, ]),
      simplify = simplify
    )
  } else {
    result <- sapply(
      1:n,
      function(i) from.theta(xx[i,]),
      simplify = simplify)
  }
  return(result)
}
