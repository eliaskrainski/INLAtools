#' A `cgeneric` model described in [cgeneric()].
setClass(
  "cgeneric",
  slots = c(f="list"),
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)
#' `rgeneric` class to define a [INLA::rgeneric()] latent model
setClass(
  "rgeneric",
  slots = c(f="list"),
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)
