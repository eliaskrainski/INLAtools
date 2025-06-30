#' Organize data for the latent GMRF C interface for `INLA`.
setClass(
  "cgeneric",
  slots = c(f = "list"),
  validity = function(object) {
    d <- b <- FALSE
    a <- all(c("model", "n", "cgeneric") %in%
               names(object$f))
    if(a) {
      b <- all(c("model", "shlib", "n", "debug", "data")
               %in% names(object$f$cgeneric))
    }
    if(a & b) {
      d <- all(c("ints", "doubles", "characters", "matrices", "smatrices")
               %in% names(object$f$cgeneric$data))
    }
    a & b & d
  }
)
#' Organize data for the latent GMRF R interface for `INLA`.
setClass(
  "rgeneric",
  slots = c(f = "list"),
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)
