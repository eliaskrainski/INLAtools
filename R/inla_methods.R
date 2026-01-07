#' Define the prec method for an inla output object
#' @param model an inla output
#' @param prior logical (default is TRUE) to indicate
#' if it should retrieve the prior or the posterior precision
#' @details
#' extract the joint prior precision for the
#' latent field at the mode of the hyperparameters
#'
#' @export
prec.inla <- function(model) {
  if(is.null(model$misc$config$config)) {
    warning("inla.rerun() with config = TRUE in control.compute.")
    model$.args$control.compute$config <- TRUE
    model <- do.call("inla", args = model$.args)
  }
  if(prior) {
    Qu <- Sparse(
      model$misc$config$config[[1]]$Qprior
    )
  } else {
    Qu <- Sparse(
      model$misc$config$config[[1]]$Q
    )
  }
  Q <- Matrix::sparseMatrix(
      i = Qu@i + 1L,
      j = Qu@j + 1L,
      x = Qu@x,
      symmetric = TRUE,
      repr = "T"
    )
  return(Q)
}

