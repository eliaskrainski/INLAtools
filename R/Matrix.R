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
#'  This is based in INLA::inla.as.sparse(),
#' but allow all combinations of 'na.rm' and 'zeros.rm'.
#' @export
Sparse <- function(A,
                   unique = TRUE,
                   na.rm = FALSE,
                   zeros.rm = FALSE) {
  if (!inherits(A, "Matrix")) {
    A <- as(A, "Matrix")
  }
  if (unique) {
    A <- as(as(as(as(A, "dMatrix"),
                  "generalMatrix"),
               "CsparseMatrix"),
            "TsparseMatrix")
  } else {
    if (!inherits(A, "dgTMatrix")) {
      A <- as(as(as(A, "dMatrix"),
                 "generalMatrix"),
              "TsparseMatrix")
    }
  }
  ## different logic than 'INLA::inla.as.sparse()'
  if (na.rm) {
    x.na <- is.na(A@x)
    if (any(x.na)) {
      idx.na <- which(x.na)
      A@x <- A@x[-idx.na]
      A@i <- A@i[-idx.na]
      A@j <- A@j[-idx.na]
    }
  }
  if (zeros.rm) {
    x.zero <- is.zero(A@x) ## changed from original
    if (any(x.zero)) {
      idx.zero <- which(x.zero)
      A@x <- A@x[-idx.zero]
      A@i <- A@i[-idx.zero]
      A@j <- A@j[-idx.zero]
    }
  }
  return(A)
}
#' Padding (a list of) sparse matrices.
#' @param M 'Matrix' (or a list of them).
#' @param relative logical. If 'M" is a list,
#' it indicates if it is to be returned a relative index
#' and the value for each matrix. See details.
#' @param ... additional arguments passed to [Sparse].
#' @returns If a unique matrix is given, return the
#' upper triangle considering the 'T' representation
#' in the `dgTMatrix`, from the `Matrix` package.
#' If a list of matrices is given,
#' return a list of two elements: 'graph' and 'xx'.
#' The 'graph' is the union of the graph from each matrix.
#' If relative=FALSE, 'xx' is a matrix with number of column equals
#' the the number of matrices imputed.
#' If relative=TRUE, it is a list of length equal the number
#' of matrices imputed. See details.
#' @details
#' This is useful to prepare a matrix, or a list of,
#' sparse matrices for use in some 'cgeneric' code.
#'
#' Define a graph of the union of the supplied matrices
#' and return the row ordered diagonal plus upper triangle
#' after padding with zeroes each one so that
#' all the returned matrices have the same pattern.
#'
#' If relative=FALSE, each columns of 'xx' is the
#' elements of the corresponding matrix after being padded
#' to fill the pattern of the union graph.
#' If relative=TRUE, each element of 'xx' would be a list
#' with a relative index, 'r', for each non-zero elements of
#' each matrix is returned relative to the union graph,
#' the non-lower elements, 'x', of the corresponding matrix,
#' and a vector, 'o', with the number of non-zero elements
#' for each line of each resulting matrix.
#' @export
#' @examples
#' A <- sparseMatrix(
#'   i = c(1, 1, 2, 3, 3, 5),
#'   j = c(2, 5, 3, 4, 5, 5),
#'   x = -c(0:3,NA,1), symmetric = TRUE)
#' A
#' upperPadding(A)
#' upperPadding(A, na.rm = TRUE)
#' upperPadding(A, zeros.rm = TRUE)
#' upperPadding(A, na.rm = TRUE, zeros.rm = TRUE)
#' B <- Diagonal(nrow(A), -colSums(A, na.rm = TRUE))
#' B
#' upperPadding(list(a = A, b = B), na.rm = TRUE, zeros.rm = TRUE)
#' upperPadding(list(a = A, b = B), relative = TRUE)
upperPadding <-
  function(M, relative = FALSE, ...) {
    .uof <- function(m) {
      return(intersect(order(m@i), which(m@j >= m@i)))
    }
    if (is(M, "list")) {
      M <- lapply(M, Sparse,  ...)
      graph <- Sparse(Reduce(
        "+",
        lapply(M, function(m) {
          m@x <- m@x * 0.0 + 1.0
          return(m)
        })
      ))
      uo <- .uof(graph)
      xx <- sapply(M, function(m) {
        m <- graph * 0 + m
        return(m@x[uo])
      })
      graph@i <- graph@i[uo]
      graph@j <- graph@j[uo]
      graph@x <- rep(1, length(uo))
      if (relative) {
        xx <- apply(xx, 2, function(x) {
          r <- which(x != 0)
          return(list(
            r = r, x = x[r]
          ))
        })
      }
      return(list(graph = graph, xx = xx))
    } else {
      M <- upperPadding(list(M))
      M$graph@x <- drop(M$xx)
      return(M$graph)
    }
  }
