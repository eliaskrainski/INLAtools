#' @describeIn methods
#' Define the prec method for an inla output object
#' @export
prec.inla <- function(model, ...) {
  if(is.null(model$misc$config$config)) {
    warning("inla.rerun() with config = TRUE in control.compute.")
    model$.args$control.compute$config <- TRUE
    model <- do.call("inla", args = model$.args)
  }
  Qu <- Sparse(
    model$misc$config$config[[1]]$Qprior
  )
  Q <- Matrix::sparseMatrix(
      i = Qu@i + 1L,
      j = Qu@j + 1L,
      x = Qu@x,
      symmetric = TRUE,
      repr = "T"
    )
  return(Q)
}

