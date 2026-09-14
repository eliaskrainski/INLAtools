#' function to make Y: for multiple likelihood with INLA
#' @param data a `data.frame` or matrix.
#' @param by.var logical to indicate if the order is by variable
#' @details
#' extract the joint prior precision for the
#' latent field at the mode of the hyperparameters
#' @examples
#' dataf <- data.frame(y1 = 1:3, y2 = rnorm(3))
#' data.expand(dataf, by.var = FALSE)
#' data.expand(dataf, by.var = TRUE)
#' @export
data.expand <- function(data, by.var = FALSE) {
  stopifnot((ni <- nrow(data))>1)
  stopifnot((nj <- ncol(data))>1)
  if(by.var) {
    out <- sapply(1:nj, function(j)
      c(rep(NA, ni*(j-1)), data[, j],
        rep(NA, ni*(nj-j))))
  } else {
    nd <- ni * nj
    out <- kronecker(as.matrix(data), matrix(1,nj,1))
    for(k in 1:nj) {
      ina <- setdiff(1:nd, seq(k, nd, nj))
      out[ina, k] <- NA
    }
  }
  colnames(out) <- colnames(data)
  return(out)
}
