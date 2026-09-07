#' Define the method to extract the precision from an inla output object.
#' @param model an inla output
#' @param ... used to pass the 'prior' argument,
#' as logical (default is TRUE) to indicate if it is
#' to retrieve the prior or the posterior precision,
#' 'id.config' to inform which config it is to be taken,
#' if not provided it will take the one with highest
#' log-posterior.
#' @details
#' extract the joint prior precision for the
#' latent field at the mode of the hyperparameters
#'
#' @export
inlaQ <- function(model, ...) {
  if(is.null(model$misc$config$config)) {
    warning("inla.rerun() with config = TRUE in control.compute.")
    model$.args$control.compute$config <- TRUE
    model <- do.call("inla", args = model$.args)
  }
  dotArgs <- list(...)
  prior <- dotArgs$prior
  if(is.null(prior))
    prior <- TRUE
  if(any(names(dotArgs)=="id.config")) {
    id.config <- dotArgs$id.config[1]
  } else {
    if(model$.args$control.inla$int.strategy=='ccd') {
      id.config <- 1
    } else {
      id.config <- which.max(
        sapply(model$misc$config$configs,
               function(g) g$log.posterior))
    }
  }
  if(prior) {
    Qu <- Sparse(
      model$misc$config$config[[id.config]]$Qprior
    )
  } else {
    Qu <- Sparse(
      model$misc$config$config[[id.config]]$Q
    )
  }
  Q <- Matrix::sparseMatrix(
      i = Qu@i + 1L,
      j = Qu@j + 1L,
      x = Qu@x,
      symmetric = TRUE
    )
  return(Sparse(Q))
}

