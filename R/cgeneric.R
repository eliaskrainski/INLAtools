#' @rdname cgeneric-class
#' @description
#' A GMRF is defined from model parameters \eqn{\theta} that
#' would parametrize a (sparse) precision matrix.
#'
#' The elements of a GMR are:
#'  *  `graph` to define the non-zero precision matrix pattern.
#'  only the upper triangle including the diagonal is needed.
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
#' @note
#' The `graph` and `Q` non-zero pattern should match,
#' its elements should be ordered by row,
#'and only its upper part stored.
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
#' @rdname cgeneric-class
#' @param model object class for what a `cgeneric` method exists.
#' E.g., if it is a character, a specific function will be called.
#' E.g. cgeneric("iid", ...") calls cgeneric_iid(...).
#' @param debug integer, used as verbose in debug.
#' @param useINLAprecomp logical, indicating if it is to use
#' the shared object previously copied and compiled by INLA.
#' @param package character giving the name of the package
#' that contains the `cgeneric` model.
#' @param libpath character, to inform the full path to the
#'  shared dynamic library object (this override the
#'  arguments `useINLAprecomp` and `package`).
#' @param ... additional arguments passed to to methods.
#' Some arguments can be used to define specific behavior,
#' such as `debug` (integer, used as verbose in debug),
#' `useINLAprecomp` (logical, indicating if it is to use
#' the shared object previously copied and compiled by INLA),
#' `package` (character used if `useINLAprecomp` is TRUE, with
#' the package name to build the path) and `libpath` (character,
#'  with the path to the shared dynamic library object: this
#'  override `useINLAprecomp` and `package`).
#' @importFrom methods is
#' @export
cgeneric.character <- function(
    model,
    debug = FALSE,
    package,
    useINLAprecomp = TRUE,
    libpath = NULL,
    ...) {

  ## some "convention" here:
  if(length(grep("cgeneric_", model))==0) {
    model <- paste0("cgeneric_", model)
  }
  if(any(model %in%
         paste0("cgeneric_",
                c("iid", "generic0", "libpath")))) {
    if(debug) {
      cat("call", model, "\n")
    }
    return(do.call(
      what = model,
      args = c(
        list(debug = debug),
        list(...)
      )
    ))
  }

  if(is.null(libpath)) {
    libpath <- cgeneric_libpath(
      fName = model,
      debug = debug,
      package = package,
      useINLAprecomp = useINLAprecomp
    )
  }

  ## do what INLA::inla.cgeneric.define() does

  d.args <- list(...)
 ## d.args <- eval(as.list(match.call())[-c(1:2)])
  nargs <- names(d.args)
  if(any(nargs == ""))
    stop("Please name the arguments!")
  if(!any(nargs == "n"))
    stop("Please provid 'n'!")
  n <- as.integer(d.args$n)
  stopifnot(n>=1)
  d.args <- d.args[setdiff(1:length(d.args),
                           which(nargs=='n'))]

  args <- c(
    list(model = model,
         n = as.integer(n),
	       debug = as.integer(debug),
	       shlib = as.character(libpath)),
    d.args)

  nM <- length(iM <- which(
    sapply(args, function(x)
      is(x, "Matrix"))))
  nm <- length(im <- which(
    sapply(args, is.matrix)))
  ni <- length(ii <- setdiff(
    which(sapply(args, is.integer)), im))
  nd <- length(id <- setdiff(
    which(sapply(args, is.double)), im))
  nc <- length(ic <- which(
    sapply(args, is.character)))
  stopifnot(ni>1)
  stopifnot(ni>1)

  if(debug) {
    cat("The cgeneric model data contains:\n",
        ni, "ints",
        nc, "characters",
        nd, "doubles",
        nm, "matrices, and",
        nM, "smatrices\n")
  }

  cmodel <- args[c("model", "shlib", "n", "debug")]
  cmodel$data <- vector("list", 5L)
  names(cmodel$data) <- c(
    "ints", "doubles", "characters",
    "matrices", "smatrices")
  cmodel$data$ints <- args[ii]
  if(nd>0) cmodel$data$doubles <- args[id]
  cmodel$data$characters <- args[ic]
  if(nm>0) {
    cmodel$data$matrices <- args[im]
    for(i in 1:nm) {
      cmodel$data$matrices[[i]] <-
        c(dim(cmodel$data$matrices[[i]]),
          cmodel$data$matrices[[i]])
    }
  }
  if(nM>0) {
	  cmodel$data$smatrices <- args[iM]
	  for(i in 1:nM) {
	    smi <- upperPadding(cmodel$data$smatrices[[i]])
	    cmodel$data$smatrices[[i]] <- c(
	      dim(smi), length(smi@x),
	      smi@i, smi@j, smi@x
	    )
	  }
  }
  class(cmodel) <- c("cgeneric",
		     "inla.cgeneric") ## this is needed in INLA::f()
  cmodel <- list(f=list(
    model = "cgeneric",
    n = cmodel$n,
    cgeneric = cmodel))
  class(cmodel) <- class(cmodel$f$cgeneric)
  return(cmodel)
}
#' @rdname cgeneric-class
#' @param fName character with the name of the `cgeneric`
#' builder function. This argumen is used by `cgeneric_libpath`
#' to build the shared lib path.
#' @returns `cgeneric_libpath' returns a character with
#' the path to the shared lib.
cgeneric_libpath <- function(
    fName,
    package,
    useINLAprecomp = FALSE,
    debug = FALSE) {

  if(missing(package)) {
    package <- attr(
      findGetFunction(
      fName = fName,
      debug = debug
    ), "package")
  }

  nbit <- 8 * .Machine$sizeof.pointer
  if(useINLAprecomp) {
    OS <- .Platform$OS.type
    if(OS=="unix") {
      OSb <- paste0("linux/", nbit, "bit/")
      if(!is.na(file.info("/Library")$isdir)) {
        OSb <- paste0("mac/", nbit, "bit/")
      }
      if(Sys.info()[["machine"]] == "arm64") {
        OSb <- "mac.arm64/"
      }
    } else {
      OSb <- paste0(OS, "/", nbit, "bit/")
    }
    shlib <- paste0(
      find.package("INLA"), "/bin/", OSb,
      "external/", package,
      "/lib", package, ".so")
    if(debug) {
      cat("INLA compiled shared lib at:\n",
          shlib, "\n")
    }
  } else {
    libpath <- system.file("libs",
                           package = package)
    if (Sys.info()["sysname"] == "Windows") {
      shlib <- file.path(
        libpath,
        paste0("x", nbit, "/", package, ".dll"))
    } else {
      shlib <- file.path(libpath,
                         paste0(package, ".so"))
    }
    if(debug) {
      cat(package, "compiled shared lib at:\n",
          shlib, "\n")
    }
  }
  return(normalizePath(shlib))
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
