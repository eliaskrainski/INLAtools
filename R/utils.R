#' To check package version and load
#' @param name character with the name of the package
#' @param minimum_version character with the minimum required version
#' @param quietly logical indicating if messages shall be printed
#' @note
#' Original in inlabru package function check_package_version_and_load
#' @export
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
#' @returns function.
#' The (first) package name where it was found
#' is returned as an attribute named "package"
#' @details
#' if 'missing(package)' it will search on the loaded
#' packages, first in the exported functions, and then
#' among the non-exported ones.
#' NOTE: 'package' can include any installed package,
#' see [installed.packages()]
#' @export
findGetFunction <- function(fName, package, debug = FALSE) {

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
