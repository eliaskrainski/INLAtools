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
#'  The original code is in inla.as.sparse()
#'  function of the 'INLA' package.
Sparse <- function(A, unique = TRUE, na.rm = FALSE, zeros.rm = FALSE) {
  if (!inherits(A, "Matrix")) {
    A <- as(A, "Matrix")
  }
  if (unique) {
    A <- as(as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
            "TsparseMatrix")
  }
  else {
    if (!inherits(A, "dgTMatrix")) {
      A <- as(as(as(A, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    }
  }
  if (na.rm) {
    idx.na <- is.na(A@x)
    if (any(idx.na)) {
      A@x[idx.na] <- 0
    }
  }
  if (zeros.rm) {
    x.zero <- (A@x == 0)
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
#' @param pkg character with the name of the package
#' @param minimum_version character with the minimum required version
#' @param quietly logical indicating if messages shall be printed
#' @note
#' The original code is in check_package_version_and_load()
#' function of the 'inlabru' package
pkgCheck <-
  function(pkg, minimum_version, quietly = FALSE) {
    version <- tryCatch(utils::packageVersion(pkg),
                        error = function(e) NA_character_
    )
    if (is.na(version)) {
      if (!quietly) {
        message(paste0("Package '", pkg, "' is not installed."))
      }
      return(NA_character_)
    }
    if (version < minimum_version) {
      if (!quietly) {
        message(paste0(
          "Installed '", pkg, "' version is ", version, " but ",
          "version >= ", minimum_version, " is required."
        ))
      }
      return(NA_character_)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!quietly) {
        message("Package '", pkg, "' not loaded safely.")
      }
      return(NA_character_)
    }
    return(version)
  }
