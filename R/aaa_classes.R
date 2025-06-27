#' The `cgeneric` class organize data needed
#' to work with Gaussian Markov Random Fields - GMRF,
#' defined as a C interface for `INLA`.
setClass(
  "cgeneric",
  slots = c(f = "list"),
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)
#' The `rgeneric` class organize data needed
#' to work with Gaussian Markov Random Fields - GMRF,
#' defined as a R interface for `INLA`.
setClass(
  "rgeneric",
  slots = c(f = "list"),
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)
