#' To store in i,j,x sparse matrix format
#' @param A matrix or Matrix
#' @param unique logical (default is TRUE) to
#' ensure that the internal representation
#' is unique and there are no duplicated entries.
#' (Do not change this unless you know what you are doing.)
#' @param na.rm logical (default is FALSE) indicating
#' if it is to replace ‘NA’'s in the matrix with zeros.
#' @param zeros.rm logical (default is FALSE)
#' indicating if it is to remove zeros in the
#' matrix. Applied after `na.rm`.
#' @note
#'  The original code is INLA::inla.as.sparse()
Sparse <- function(A,
                   unique = TRUE,
                   na.rm = FALSE,
                   zeros.rm = FALSE) {
  if (!inherits(A, "Matrix")) {
    A <- as(A, "Matrix")
  }
  if (unique) {
    A <- as(as(as(as(A, "dMatrix"),
                  "generalMatrix"),
               "CsparseMatrix"),
            "TsparseMatrix")
  } else {
    if (!inherits(A, "dgTMatrix")) {
      A <- as(as(as(A, "dMatrix"),
                 "generalMatrix"),
              "TsparseMatrix")
    }
  }
  if (na.rm) {
    idx.na <- is.na(A@x)
    if (any(idx.na)) {
      A@x[idx.na] <- 0
    }
  }
  if (zeros.rm) {
    x.zero <- is.zero(A@x) ## changed from original
    if (any(x.zero)) {
      idx.zero <- which(x.zero)
      A@x <- A@x[-idx.zero]
      A@i <- A@i[-idx.zero]
      A@j <- A@j[-idx.zero]
    }
  }
  return(A)
}
#' To check package version and load
#' @param name character with the name of the package
#' @param minimum_version character with the minimum required version
#' @param quietly logical indicating if messages shall be printed
#' @note
#' The original code is inlabru:::check_package_version_and_load()
packageCheck <- function(name, minimum_version, quietly = FALSE) {
  version <- tryCatch(utils::packageVersion(name),
                      error = function(e) NA_character_
  )
  if (is.na(version)) {
    if (!quietly) {
      message(paste0("Package '", name, "' is not installed."))
    }
    return(NA_character_)
  }
  if (version < minimum_version) {
    if (!quietly) {
      message(paste0(
        "Installed '", name, "' version is ", version, " but ",
        "version >= ", minimum_version, " is required."
      ))
    }
    return(NA_character_)
  }
  if (!requireNamespace(name, quietly = TRUE)) {
    if (!quietly) {
      message("Package '", name, "' not loaded safely.")
    }
    return(NA_character_)
  }
  return(version)
}
#' Search a function and retrieve it.
#' @param fName character with the name of the function
#' @param package character with the package name
#' @param debug logical indicating if it is to print
#' intermediate progress finding
#' @returns function
#' @details
#' The (first) package name where it was found
#' is returned as an atribute named "package"
#' @export
findAndGetFunction <- function(fName, package, debug = FALSE) {

  if(missing(package)) {
    pkgs <- .packages()
  } else {
    pkgs <- package
    package <- NULL
  }

  ## first try exported
  for(i in 1:length(pkgs)) {
    afn <- paste0(pkgs[i], "::", fName)
    efn <- try(eval(str2lang(afn)), silent = TRUE)
    if(is.function(efn)) {
      package <- pkgs[i]
      if(debug) {
        cat("Found", fName, "in", package, "!\n")
      }
      break
    }
  }

  if(is.null(package)) {
    ## try non-exported
    for(i in 1:length(pkgs)) {
      afn <- paste0(pkgs[i], ":::", fName)
      efn <- try(eval(str2lang(afn)), silent = TRUE)
      if(is.function(efn)) {
        package <- pkgs[i]
        if(debug) {
          cat(fName, "non-exported from", package, "!\n")
        }
        break
      }
    }
  }

  if(is.null(package)) {
    stop("Function not found in the tried packages!")
  }

  attr(efn, "package") <- package
  return(efn)
}
