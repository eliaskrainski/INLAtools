#' The `cgeneric` class for [cgeneric()].
setClass(
  "cgeneric",
  slots = c(f = "list"),
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)
#' The `rgeneric` class for [rgeneric()].
setClass(
  "rgeneric",
  slots = c(f = "list"),
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)
