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
#' @returns a method to build a `cgeneric` should return
#' a named list of `cgeneric` class that contains a
#'  named list `f` that contains (at least):
#'  * `model` a character always equal to `cgeneric`,
#'  * `n` an integer greater than 0, and
#'  * `cgeneric` as a named list that contains the
#'  data needed to define the model. Each element on
#'  ...$f$cgeneric is also a named list containing
#'  `ints`, `doubles`, `characters`, `matrices`
#'  and `smatrices`.
#'  * (possible) `extraconstr` as a named list with: `A` as a
#'   `n` times `k` matrix and `e` as a length `k` vector.
#'  * (possible) `bm_mapper` (TO DO) mapper for `inlabru` package.
#'
#'  The `cgeneric_shlib` function returns a `character`
#'  with the path to the shared lib.
#' @seealso [INLA::cgeneric()] and [INLAtools-methods()]
#' @export
#' @example demo/cgeneric.R
cgeneric <- function(model, ...) {
  UseMethod("cgeneric")
}
#' @rdname cgeneric-class
#' @param model object class for what a `cgeneric` method exists.
#' E.g., if it is a character, a specific function will be called.
#' E.g. cgeneric("iid", ...") calls cgeneric_iid(...).
#' @param ... arguments passed from the
#' [cgeneric()] methods. FIt should include
#' `n` and `debug`. For `cgenericBuild` it should
#'  `model` as a character string with the name of the
#'  C function and `shlib` as the path to the
#'  shared object containing such function.
#'  If `shlib` is not provided it can be built
#'  using `inla_shlib` from the arguments,
#'  `package` (character with the R package containing it),
#'  `useINLAprecomp` (logical to indicate if `INLA` contains it
#'  and to use it).
#' @importFrom methods is
#' @export
cgeneric.character <- function(
    model,
    ...) {

  dotArgs <- list(...)

  ## some "convention" here:
  if(length(grep("cgeneric_", model))==0) {
    model <- paste0("cgeneric_", model)
  }
  if(any(model %in%
         paste0("cgeneric_",
                c("iid", "generic0", "shlib")))) {
    if(!is.null(dotArgs$debug) && dotArgs$debug) {
      cat("call", model, "\n")
    }
    return(do.call(
      what = model,
      args = dotArgs
      )
    )
  } else {

    out <- try(do.call(
      what = model,
      args = dotArgs
    ), silent = TRUE)

    if((inherits(out, "cgeneric")) |
       (inherits(out, "inla.cgeneric"))) {
      return(out)
    } else {

      fn <- findGetFunction(
        fName = model
      )
      if(is.function(fn)) {
        return(cgeneric(model = fn, ...))
      }
    }
  }

}
#' @rdname cgeneric-class
#' @export
cgeneric.function <- function(
    model,
    ...) {
  out <- do.call(
    what = model,
    args = list(...)
    )
  return(out)
}


#' @describeIn cgeneric-class Returns the model object unchanged.
#' @export
cgeneric.cgeneric <- function(model, ...) {
  return(model)
}

#' @describeIn cgeneric-class Converts a regular `inla.cgeneric` object to `cgeneric`.
#' @export
cgeneric.inla.cgeneric <- function(model, ...) {
  # TODO: Is it enough to just add the "cgeneric" class name, or does the object need to be
  # regenerated/modified as well?
  warning("TODO: check whether basic inla.cgeneric objects fulfil the assumptions of the cgeneric class")
  class(model) <- c("cgeneric", class(model))
  return(model)
}


#' @rdname cgeneric-class
#' @export
cgenericBuilder <- function(
    ...) {

  ## do what INLA::inla.cgeneric.define() does

  dotArgs <- list(...)
  nargs <- names(dotArgs)
  if(any(nargs == ""))
    stop("Please name the arguments!")
  if(!any(nargs == "n"))
    stop("Please provid 'n'!")
  if(!any(nargs == "debug"))
    dotArgs$debug <- FALSE
  if(!any(nargs == "shlib"))
    stop("Please provid 'shlib'!")
  if(!any(nargs == "model"))
    stop("Please provid 'model'!")
  stopifnot(dotArgs$n>=1)
  ind <- pmatch(c("n", "debug"),
                names(dotArgs))
  dotArgs <- c(list(
    n = dotArgs$n,
    debug = dotArgs$debug),
    dotArgs[setdiff(1:length(dotArgs),ind)]
  )
  dotArgs$n <- as.integer(dotArgs$n)
  dotArgs$debug <- as.integer(dotArgs$debug)
  dotArgs$model <- as.character(dotArgs$model)
  dotArgs$shlib <- as.character(dotArgs$shlib)

  nM <- length(iM <- which(
    sapply(dotArgs, function(x)
      is(x, "Matrix"))))
  nm <- length(im <- which(
    sapply(dotArgs, is.matrix)))
  ni <- length(ii <- setdiff(
    which(sapply(dotArgs, is.integer)), im))
  nd <- length(id <- setdiff(
    which(sapply(dotArgs, is.double)), im))
  nc <- length(ic <- which(
    sapply(dotArgs, is.character)))
  stopifnot(ni>1)
  stopifnot(ni>1)

  if(dotArgs$debug) {
    cat("The cgeneric model data contains:\n",
        ni, "ints",
        nc, "characters",
        nd, "doubles",
        nm, "matrices, and",
        nM, "smatrices\n")
  }

  cmodel <- dotArgs[c("model", "shlib", "n", "debug")]
  cmodel$data <- vector("list", 5L)
  names(cmodel$data) <- c(
    "ints", "doubles", "characters",
    "matrices", "smatrices")
  cmodel$data$ints <- dotArgs[ii]
  if(nd>0) cmodel$data$doubles <- dotArgs[id]
  cmodel$data$characters <- dotArgs[ic]
  if(nm>0) {
    cmodel$data$matrices <- dotArgs[im]
    for(i in 1:nm) {
      cmodel$data$matrices[[i]] <-
        c(dim(cmodel$data$matrices[[i]]),
          t(cmodel$data$matrices[[i]])) ## columnwise...
    }
  }
  if(nM>0) {
	  cmodel$data$smatrices <- dotArgs[iM]
	  for(i in 1:nM) {
	    smi <- upperPadding(cmodel$data$smatrices[[i]])
	    cmodel$data$smatrices[[i]] <- c(
	      dim(smi), length(smi@x),
	      smi@i, smi@j, smi@x
	    )
	  }
  }
  class(cmodel) <- c("inla.cgeneric.f", ## this is needed in INLA::f() from August 2025
                     "inla.cgeneric") ##  this is needed for older INLA
  cmodel_wrapper <- structure(
    list(
      f = list(
        model = "cgeneric",
        n = as.integer(cmodel$n),
        cgeneric = cmodel
      ),
      mapper = inlabru::bm_index(cmodel$n)
    ),
    class = c("cgeneric", "inla.cgeneric")
  )
  return(cmodel_wrapper)
}
#' @rdname cgeneric-class
#' @param debug integer, used as verbose in debug.
#' @param useINLAprecomp logical, indicating if it is to use
#' the shared object previously copied and compiled by INLA.
#' @param package character giving the name of the package
#' that contains the `cgeneric` model.
#' @export
cgeneric_shlib <- function(
    debug,
    package,
    useINLAprecomp) {

  if(missing(package) || is.null(package)) {
    stop("please provide package!")
  }
  if(missing(debug)) {
    debug <- FALSE
  }
  if(missing(useINLAprecomp)) {
    useINLAprecomp <- TRUE
  }

  nbit <- 8 * (.Machine$sizeof.pointer)
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
    shlib <- system.file("libs",
                         package = package)
    if (Sys.info()["sysname"] == "Windows") {
      shlib <- file.path(
        shlib,
        paste0("x", nbit, "/", package, ".dll"))
    } else {
      shlib <- file.path(shlib,
                         paste0(package, ".so"))
    }
    if(debug) {
      cat(package, "compiled shared lib at:\n",
          shlib, "\n")
    }
  }
  return(normalizePath(shlib))
}
#' @describeIn cgeneric-class
#' Print the cgeneric object
#' @param x a cgeneric object
#' @param ... not used
#' @export
print.cgeneric <- function(x, ...) {
  cat("cgeneric: ", x$f$cgeneric$model, ", n = ",
      x$f$cgeneric$n, ", shlib:\n", sep = "")
  cat(x$f$cgeneric$data$characters$shlib, "\n")
  d0 <- c(2L, 0L, 2L, 0, 0)
  dn <- sapply(x$f$cgeneric$data, length)
  dn0 <- dn-d0
  ii <- which(dn0>0)
  if(length(ii)>0) {
    nd <- names(dn)
    for(i in ii){
      d <- x$f$cgeneric$data[[i]]
      ndi <- names(d)
      jj <- (d0[i]+1):dn[i]
      cat(nd[i], " (", dn[i], ") : ", sep = "")
      if(i %in% c(1,2,3)) {
        cat(paste(ndi[jj], "(", sapply(d[jj], length),
                  ") ", sep = ""), "\n")
      } else {
        if(i==4) {
          dj <- sapply(d[jj], function(x) x[1:2])
          cat(paste(ndi[jj], "(", dj[1,], ",", dj[2, ],
                    ") ", sep = ""), "\n")
        }
        if(i==5) {
          dj <- sapply(d[jj], function(x) x[1:3])
          cat(paste(ndi[jj], "(", dj[1,], ",",
                    dj[2, ], ",", dj[3,],
                    ") ", sep = ""), "\n")
        }
      }
    }
  }
}
#' @describeIn cgeneric-class
#' A summary for a cgeneric object
#' @param object a cgeneric object
#' @param ... not used
#' @export
summary.cgeneric <- function(object, ...) {
  g <- graph(object)
  cat("n = ", object$f$cgeneric$n, ", graph with",
      length(g@x), "non-zeros\n", sep = "")
}
#' @describeIn cgeneric-class
#' A plot for a cgeneric object
#' @param y not used
#' @importFrom graphics image
#' @export
plot.cgeneric <- function(x, y, ...) {
  g <- graph(x)
  image(g)
}
