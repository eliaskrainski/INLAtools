#' Kronecker (product) between `cgeneric`/`rgeneric` models,
#' implemented as [kronecker()] methods.
#' @name kronecker
#' @param X `cgeneric` or `rgeneric`
#' @param Y `cgeneric` or `rgeneric`
#' @param FUN see [kronecker()]
#' @param make.dimnames see [kronecker()]
#' @param ... see [kronecker()]
#' @return if 'X' and 'Y' are `cgeneric`
#' return a `cgeneric`, else a `rgeneric`.
#' @importFrom inlabru bru_get_mapper
NULL
#> NULL

#' @rdname kronecker
#' @useDynLib INLAtools
#' @importFrom utils str
#' @export
#' @examples
#' R <- Matrix(crossprod(diff(diag(4))))
#' m1 <- cgeneric("generic0", R = R, param = c(1, NA),
#'   scale = FALSE, useINLAprecomp = FALSE)
#' m2 <- cgeneric("iid", n = 3, param = c(1, 0.5),
#'   useINLAprecomp = FALSE)
#' k21 <- kronecker(m2, m1, useINLAprecomp = FALSE)
#' prec(k21, theta = 0.0)
setMethod(
  "kronecker",
  c(X="cgeneric", Y = "cgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    if(is.null(mcall$useINLAprecomp)) {
      useINLAprecomp = TRUE
    } else {
      useINLAprecomp = mcall$useINLAprecomp
    }

    INLAvcheck <- packageCheck("INLA", "25-10-28")
    if(is.na(INLAvcheck) & useINLAprecomp) {
      useINLAprecomp <- FALSE
      warning("INLA version is old. Setting 'useINLAprecomp = FALSE'!")
    }
    shlib <- cgeneric_shlib(
      package = "INLAtools",
      useINLAprecomp = useINLAprecomp,
      debug = debug)

    cmodel <- "inla_cgeneric_kronecker"

    n1 <- as.integer(X$f$n)
    n2 <- as.integer(Y$f$n)
    N <- as.integer(n1 * n2)

    if(debug) {
      cat('n1:', n1, "n2:", n2, "n:", N, "")
    }

    ## Q1 graph index: i,j
    ij1 <- cgeneric_get(
      model = X,
      cmd = "graph",
      optimize = TRUE
    )
    stopifnot(all(ij1[[1]]<=ij1[[2]]))
    M1 <- length(ij1[[1]])
    if(debug) {
      cat('M1:', M1, "\n")
      print(str(ij1))
    }
    names(ij1) <- c("i", "j")

    ## Q2 graph index: i,j
    ij2 <- cgeneric_get(
      model = Y,
      cmd = "graph",
      optimize = TRUE
    )
    stopifnot(all(ij2[[1]]<=ij2[[2]]))
    M2 <- length(ij2[[1]])
    if(debug) {
      cat('M2:', M2, "\n")
      print(str(ij2))
    }
    names(ij2) <- c('i', 'j')

    ## For Q = Q1 (x) Q2
    ## The total number of elements in Q is
    ##  [(M1-n1)*2+n1]*[(M2-n2)*2+n2]      =
    ##    (2*M1 - n1) * (2*M2 - n2)        =
    ##  n1*n2 + 4*M1*M2 -2*M1*n2 -2*M2*n1  =
    ##  n1*n2 + 2*M1*(2*M2 - n2) - 2*M2*n1
    ## but at the upper part is made up of
    ## diagonal + half_of_the_off_diagonal
    ##  n1*n2   + M1*(2*M2 - n2) - M2*n1         =
    ##  n1*n2   + 2*M1*M2 - M1*n2 - M2*n1        =
    ##  n1*n2   + M1*M2 + M1*M2 - M1*n2 - M2*n1  =
    ##  n1*n2   + M1*M2 + M1*(M2-n2) - M2*n1     =
    ##  n1*n2   + M1*(M2-n2) + M2*(M1-n1)

    ## define u1 = M1-n1, u2 = M2-n2
    ##  n1*n2 + M1*u2        + M2*u1            =
    ##  n1*n2 + (n1 + u1)*u2 + (n2 + u2)*u1     =
    ##  n1*n2 + n1*u2 + u1*u2 + n2*u1 + u2*u1   =
    ##  n1*(n2 + u2)  + u1*(n2 + u2)  + u2*u1   =
    ##  n1*M2         + u1*M2         + u2*u1   =
    ##  (n1 + u1) * M2                + u2*u1   =
    ##         M1 * M2 + u2 * u1
    ## the i<j part from Q1
    idx1u <- which(ij1$i < ij1$j)
    u1 <- length(idx1u)
    stopifnot(u1 == (M1-n1))
    ## the i<j part from Q2
    idx2u <- which(ij2$i < ij2$j)
    u2 <- length(idx2u)
    stopifnot(u2 == (M2-n2))

    ## the number of non-diagonal elements
    M <- M1 * M2 + u1 * u2
    if(debug) {
      cat("u1:", u1, "u2:", "u2", u2, "M:", M, "\n")
    }

    ## resulting graph index i
    ii0 <- rep(ij1$i * n2, each = M2)  + ij2$i
    ## resulting graph index j
    jj0 <- rep(ij1$j * n2, each = M2)  + ij2$j
    if((u1*u2)>0) {
      ii <- c(
        ii0,
        rep(ij1$i[idx1u]*n2, each=u2) + ij2$j[idx2u])
      jj <- c(
        jj0,
        rep(ij1$j[idx1u]*n2, each=u2) + ij2$i[idx2u]
      )
    } else {
      ii <- ii0
      jj <- jj0
    }

    ## check
    stopifnot(all(ii <= jj))

    ## the order of the output should use this
    jjord <- order(jj)
    iiord <- order(ii[jjord])
    ije <- list(
      ii = ii[jjord][iiord],
      jj = jj[jjord][iiord],
      ord = jjord[iiord]
    )

    ## initial data
    ret <- structure(
      list(
        f = list(
          model = "cgeneric",
          n = as.integer(N),
          cgeneric = structure(
            list(
              model = cmodel,
              shlib = shlib,
              n = as.integer(N),
              debug = as.integer(debug)
            ),
            # inla.cgeneric is needed to support INLA before August 2025
            class = c("inla.cgeneric.f", "inla.cgeneric")
          )
        )
      ),
      class = c("cgeneric", "inla.cgeneric")
    )
    ret$f$cgeneric$data <- vector("list", 5L)
    names(ret$f$cgeneric$data) <- c(
      "ints", "doubles", "characters",
      "matrices", "smatrices")

    ## data size for each model
    ndata1 <- sapply(
      X$f$cgeneric$data,
      length)
    ndata2 <- sapply(
      Y$f$cgeneric$data,
      length)

    ### n and the data size info
      Ndata <- list(
        n = as.integer(
          c(n1, ndata1, M1,
                ndata2, M2, N, M)
        )
      )
      if(debug) {
        cat(Ndata$n, "\n")
      }
      ret$f$cgeneric$data$ints <-
        c(
          Ndata,
          ### concatenate ints from each model
          X$f$cgeneric$data$ints[-1], ## n for model 1 is already there
          Y$f$cgeneric$data$ints,
          list(
            idx1u = as.integer(idx1u - 1),
            idx2u = as.integer(idx2u - 1)
          )
        )

      if((ndata1[2]>0) | (ndata2[2]>0)) {
        ret$f$cgeneric$data$doubles <-
          c(
            X$f$cgeneric$data$doubles,
            Y$f$cgeneric$data$doubles
          )
      }

      ret$f$cgeneric$data$characters <-
      c(
        list(
          model = cmodel,
          shlib = shlib
        ),
        X$f$cgeneric$data$characters,
        Y$f$cgeneric$data$characters
      )

      if((ndata1[4]>0) | (ndata2[4]>0)) {
        ret$f$cgeneric$data$matrices <-
          c(
            X$f$cgeneric$data$matrices,
            Y$f$cgeneric$data$matrices
          )
      }

    ret$f$cgeneric$data$smatrices <-
      c(
        X$f$cgeneric$data$smatrices,
        Y$f$cgeneric$data$smatrices,
        list(Kgraph = c(
          N, N, M,
          ije$ii,
          ije$jj,
          as.numeric(ije$ord-1)
          )
        )
      )

    if((!is.null(X$f$extraconstr)) |
       (!is.null(Y$f$extraconstr))) {
        ret$f$extraconstr <-
          kronecker_extraconstr(
            X$f$extraconstr,
            Y$f$extraconstr,
            X$f$n, Y$f$n
          )
        if(debug) {
          cat(nrow(ret$f$extraconstr),
              " 'extraconstr' built!\n")
        }
    }

    # Note: kron(X,Y) gives X-major ordering (the X-index varies slowly, the
    # Y-index varies quickly), which requires the mappers to be in reverse
    # order, multi(Y,X):
    inlabruCheck <- packageCheck("inlabru", "2.13.0.9005")
    if(is.na(inlabruCheck)) {
      warning("Please install a inlabru recent version from git.")
    } else {
      ret$mapper <-
      inlabru::bm_multi(
        list(
          inlabru::bru_get_mapper(Y),
          inlabru::bru_get_mapper(X)
        )
      )
    }

    return(ret)
  }
)
#' @rdname kronecker
#' @useDynLib INLAtools
#' @importFrom utils str
#' @export
setMethod(
  "kronecker",
  c(X="cgeneric", Y = "rgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    n <- X$f$n * Y$f$n

    ini1 <- initial(X)
    nth1 <- length(ini1)
    ini2 <- initial(Y)
    nth2 <- length(ini2)

    kmodel <- function(cmd = c("graph", "Q", "mu",
                               "initial", "log.norm.const",
                               "log.prior", "quit"),
                       theta = NULL) {

      graph <- function(n, theta) {
        g1 <- cgeneric_get(X, "graph", optimize = FALSE)
        g2 <- INLA::inla.rgeneric.q(Y, "graph")
        return(kronecker(g1, g2))
      }

      Q <- function(n, theta) {
        Q1 <- cgeneric_get(X, "Q", theta = theta[1:nth1], optimize = FALSE)
        Q2 <- INLA::inla.rgeneric.q(Y, "Q", theta = theta[nth1+1:nth2], optimize = FALSE)
        QQ <- INLA::inla.as.sparse(kronecker(Q1, Q2))
        idx <- which(QQ@i <= QQ@j)
        return(QQ@x[idx])
      }

      mu <- function(n, theta)
        return(numeric(0))

      log.norm.const <- function(n, theta)
        return(numeric(0))

      log.prior <- function(n, theta) {
        return(
          prior(X, theta = theta[1:nth1]) +
            prior(Y, theta = theta[nth1+1:nth2])
        )
      }

      initial <- function(n, theta) {
        return(
          c(ini1, ini2)
        )
      }

      quit <- function(n, theta) {
        return(invisible())
      }

      cmd <- match.arg(cmd)

      ret <- do.call(
        cmd,
        args = list(n = n,
                    theta = theta
        )
      )

      return(ret)

    }

    ### follows INLA:::inla.rgeneric.define() but no assign env
    rmodel <- structure(
      list(
        f = list(
          model = "rgeneric",
          n = n,
          rgeneric = structure(
            list(
              definition =
                compiler::cmpfun(
                  kmodel,
                  options = list(optimize = 3L)),
              debug = debug,
              optimize = TRUE
            ),
            # inla.rgeneric is needed to support INLA before August 2025
            class = c("inla.rgeneric.f", "inla.rgeneric")
          )
        )
      ),
      class = c("rgeneric", "inla.rgeneric")
    )

    if((!is.null(X$f$extraconstr)) |
       (!is.null(Y$f$extraconstr))) {
      rmodel$f$extraconstr <-
        kronecker_extraconstr(
          X$f$extraconstr,
          Y$f$extraconstr,
          X$f$n, Y$f$n
        )
      if(debug) {
        cat(nrow(rmodel$f$extraconstr),
            " 'extraconstr' built\n")
      }
    }

    # Note: kron(X,Y) gives X-major ordering (the X-index varies slowly, the
    # Y-index varies quickly), which requires the mappers to be in reverse
    # order, multi(Y,X):
    inlabruCheck <- packageCheck("inlabru", "2.13.0.9005")
    if(is.na(inlabruCheck)) {
      warning("Please install a inlabru recent version from git.")
    } else {
      rmodel$mapper <-
      inlabru::bm_multi(
        list(
          inlabru::bru_get_mapper(Y),
          inlabru::bru_get_mapper(X)
        )
      )
    }

    return(rmodel)
  }
)
#' @rdname kronecker
#' @useDynLib INLAtools
#' @importFrom utils str
#' @export
setMethod(
  "kronecker",
  c(X="rgeneric", Y = "cgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    n <- X$f$n * Y$f$n

    ini1 <- initial(X)
    nth1 <- length(ini1)
    ini2 <- initial(Y)
    nth2 <- length(ini2)

    kmodel <- function(cmd = c("graph", "Q", "mu",
                               "initial", "log.norm.const",
                               "log.prior", "quit"),
                       theta = NULL) {

      graph <- function(n, theta) {
        g1 <- INLA::inla.rgeneric.q(rmodel = X, cmd = "graph")
        g2 <- cgeneric_get(Y, "graph", optimize = FALSE)
        return(kronecker(g1, g2))
      }

      Q <- function(n, theta) {
        Q1 <- INLA::inla.rgeneric.q(rmodel = X, cmd = "Q", theta = theta[1:nth1])
        Q2 <- cgeneric_get(Y, "Q", theta = theta[nth1+1:nth2], optimize = FALSE)
        QQ <- INLA::inla.as.sparse(kronecker(Q1, Q2))
        idx <- which(QQ@i <= QQ@j)
        return(QQ@x[idx])
      }

      mu <- function(n, theta)
        return(numeric(0))

      log.norm.const <- function(n, theta)
        return(numeric(0))

      log.prior <- function(n, theta) {
        return(
          INLA::inla.rgeneric.q(rmodel = X, cmd = "log.prior", theta = theta[1:nth1]) +
          prior(Y, theta = theta[nth1+1:nth2])
        )
      }

      initial <- function(n, theta) {
        return(c(ini1, ini2))
      }

      quit <- function(n, theta) {
        return(invisible())
      }

      cmd <- match.arg(cmd)

      ret <- do.call(
        cmd,
        args = list(n = n,
                    theta = theta
        )
      )

      return(ret)

    }

    ### follows INLA:::inla.rgeneric.define() but no assign env
    rmodel <- structure(
      list(
        f = list(
          model = "rgeneric",
          n = n,
          rgeneric = structure(
            list(
              definition =
                compiler::cmpfun(
                  kmodel,
                  options = list(optimize = 3L)),
              debug = debug,
              optimize = TRUE
            ),
            # inla.rgeneric is needed to support INLA before August 2025
            class = c("inla.rgeneric.f", "inla.rgeneric")
          )
        )
      ),
      class = c("rgeneric", "inla.rgeneric")
    )

    if((!is.null(X$f$extraconstr)) |
       (!is.null(Y$f$extraconstr))) {
      rmodel$f$extraconstr <-
        kronecker_extraconstr(
          X$f$extraconstr,
          Y$f$extraconstr,
          X$f$n, Y$f$n
        )
      if(debug) {
        cat(nrow(rmodel$f$extraconstr),
            " 'extraconstr' built!\n")
      }
    }

    # Note: kron(X,Y) gives X-major ordering (the X-index varies slowly, the
    # Y-index varies quickly), which requires the mappers to be in reverse
    # order, multi(Y,X):
    inlabruCheck <- packageCheck("inlabru", "2.13.0.9005")
    if(is.na(inlabruCheck)) {
      warning("Please install a inlabru recent version from git.")
    } else {
      rmodel$mapper <-
      inlabru::bm_multi(
        list(
          inlabru::bru_get_mapper(Y),
          inlabru::bru_get_mapper(X)
        )
      )
    }

    return(rmodel)
  }
)
#' @rdname kronecker
#' @importFrom utils str
setMethod(
  "kronecker",
  c(X="rgeneric", Y = "rgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    n <- X$f$n * Y$f$n

    ini1 <- initial(X)
    nth1 <- length(ini1)
    ini2 <- initial(Y)
    nth2 <- length(ini2)

    kmodel <- function(cmd = c("graph", "Q", "mu",
                               "initial", "log.norm.const",
                               "log.prior", "quit"),
                       theta = NULL) {

      graph <- function(n, theta) {
        g1 <- INLA::inla.rgeneric.q(X, "graph")
        g2 <- INLA::inla.rgeneric.q(Y, "graph")
        return(kronecker(g1, g2))
      }

      Q <- function(n, theta) {
        Q1 <- INLA::inla.rgeneric.q(rmodel = X, cmd = "Q", theta = theta[1:nth1])
        Q2 <- INLA::inla.rgeneric.q(rmodel = Y, cmd = "Q", theta = theta[nth1+1:nth2])
        QQ <- INLA::inla.as.sparse(
          kronecker(Q1, Q2))
        idx <- which(QQ@i <= QQ@j)
        return(QQ@x[idx])
      }

      mu <- function(n, theta)
        return(numeric(0))

      log.norm.const <- function(n, theta)
        return(numeric(0))

      log.prior <- function(n, theta) {
        lp1 <- INLA::inla.rgeneric.q(rmodel = X, cmd = "log.prior", theta = theta[1:nth1])
        lp2 <- INLA::inla.rgeneric.q(rmodel = Y, cmd = "log.prior", theta = theta[nth1+1:nth2])
        return(lp1 + lp2)
      }

      initial <- function(n, theta) {
        return(c(ini1, ini2))
      }

      quit <- function(n, theta) {
        return(invisible())
      }

      cmd <- match.arg(cmd)

      ret <- do.call(
        cmd,
        args = list(
          n = n,
          theta = theta
        )
      )
      return(ret)
    }

    rmodel <- INLA::inla.rgeneric.define(
      model = kmodel,
      optimize = TRUE
    )
    class(rmodel) <- c("rgeneric", class(rmodel))

    if((!is.null(X$f$extraconstr)) |
       (!is.null(Y$f$extraconstr))) {
      rmodel$f$extraconstr <-
        kronecker_extraconstr(
          X$f$extraconstr,
          Y$f$extraconstr,
          X$f$n, Y$f$n
        )
      if(debug) {
        cat(nrow(rmodel$f$extraconstr),
            " 'extraconstr' built!\n")
      }
    }

    # Note: kron(X,Y) gives X-major ordering (the X-index varies slowly, the
    # Y-index varies quickly), which requires the mappers to be in reverse
    # order, multi(Y,X):
    inlabruCheck <- packageCheck("inlabru", "2.13.0.9005")
    if(is.na(inlabruCheck)) {
      warning("Please install a inlabru recent version from git.")
    } else {
      rmodel$mapper <-
      inlabru::bm_multi(
        list(
          inlabru::bru_get_mapper(Y),
          inlabru::bru_get_mapper(X)
        )
      )
    }

    return(rmodel)
  }
)
#' Kronecker (product) between `extraconstr`,
#' implemented for [kronecker()] methods.
#' @name extraconstr
#' @param c1,c2 named list with two elements:
#' `A` and `e`, where `nrow(A)` should be equal
#' to `length(e)`. These are constraint definitions.
#' @param n1,n2 integer with each model's length.
#' @returns The constraint definition for the
#' whole latent model built from the Kronecker product.
#' A length two named list. 'A' a matrix and
#' 'e' a vector where nrow(A)=length(e) and
#' ncol(A)=(n1*n2).
kronecker_extraconstr <- function(c1, c2, n1, n2) {
  if(is.null(c1)) {
    if(is.null(c2)) {
      return(NULL)
    } else {
      ret <- list(
        A = kronecker(diag(n1), c2$A),
        e = rep(c2$e, n1)
      )
    }
  } else {
    if(is.null(c2)) {
      ret <- list(
        A = kronecker(c1$A, diag(n2)),
        e = rep(c1$e, each = n2)
      )
    } else {
      ret <- list(
        A = rbind(
          kronecker(c1$A, diag(ncol(c2$A))),
          kronecker(diag(ncol(c1$A)), c2$A)
        )[-1, , drop = FALSE], ## remove one redundant
        e = c(rep(c1$e, each = ncol(c2$A)),
              rep(c2$e, ncol(c1$A)))[-1] ## rm 1
      )
    }
  }
  return(ret)
}


#' @title Combine two or more `cgeneric` or `rgeneric` models
#' @description Constructs a multiple kronecker product model from a list of model
#' objects. The resulting model contains a corresponding [inlabru::bm_multi()]
#' mapper. This can be used as an alternative to a binary tree of kronecker
#' product models.
#' @param models A list of `cgeneric` or `rgeneric` models, optionally with names
#' @details The last model in the list has the slowest index variation, and the
#'   first model has the fastest index variation. This matches the latent
#'   variable ordering of standard `INLA:f()` model components with
#'   `(main, group, replicate)`.
#' @param \dots Arguments passed on to every `kronecker()` call.
#' @returns A 'cgeneric' or 'rgeneric' model object, containing a
#'   multi-kronecker product model, with a corresponding [inlabru::bm_multi()]
#'   mapper.
#' @rdname multi_generic_model
#' @export
#' @examples
#' R1 <- Matrix(crossprod(diff(diag(4))))
#' m1 <- cgeneric("generic0", R = R1, param = c(1, NA),
#'   scale = FALSE, useINLAprecomp = FALSE)
#' R2 <- Matrix(crossprod(diff(diag(3))))
#' m2 <- cgeneric("generic0", R = R2, param = c(1, NA),
#'   scale = FALSE, useINLAprecomp = FALSE)
#' m3 <- cgeneric("iid", n = 2, param = c(1, 0.5),
#'   useINLAprecomp = FALSE)
#' multi123 <- multi_generic_model(
#'   list(m1 = m1, m2 = m2, m3 = m3),
#'   useINLAprecomp = FALSE
#' )
#' prec(multi123, theta = 0.0)
#' if(!is.na(packageCheck("inlabru", "2.13.0.9005"))) {
#'   print(multi123$mapper)
#' }
multi_generic_model <- function(models, ...) {
  stopifnot(is.list(models))
  stopifnot(length(models) >= 1L)

  models <- lapply(models, function(model) {
    if (inherits(model, c("rgeneric", "inla.rgeneric"))) {
      rgeneric(model)
    } else if (inherits(model, c("cgeneric", "inla.cgeneric"))) {
      cgeneric(model)
    } else {
      stop("Each 'models' element must be convertible to class 'cgeneric' or 'rgeneric'")
    }
  })

  ret <- models[[1]]
  for (i in seq_along(models)[-1]) {
    ret <- Matrix::kronecker(models[[i]], ret, ...)
  }

  inlabruCheck <- packageCheck("inlabru", "2.13.0.9005")
  if(is.na(inlabruCheck)) {
    warning("Please install a inlabru recent version from git.")
  } else {
    ret[["mapper"]] <- inlabru::bm_multi(
      lapply(models, inlabru::bru_get_mapper)
    )
  }

  return(ret)
}
