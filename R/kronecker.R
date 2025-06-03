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
    if(is.null(mcall$libpath)) {
      libpath <- NULL
    } else {
      libpath <- mcall$libpath
    }

    model <- "inla_cgeneric_kronecker"
    if (is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("INLAtools")
      } else {
        libpath <- system.file("libs", package = "INLAtools")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "x64/INLAtools.dll")
        } else {
          libpath <- file.path(libpath, "INLAtools.so")
        }
      }
    }

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
    ret <- list(
      f = list(
        model = "cgeneric",
        n = as.integer(N),
        cgeneric = list(
          model = model,
          shlib = libpath,
          n = as.integer(N),
          debug = as.integer(debug),
          data = list()
        )
      )
    )

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

    ret$f$cgeneric$data$doubles <-
      c(
        X$f$cgeneric$data$doubles,
        Y$f$cgeneric$data$doubles
      )

    ret$f$cgeneric$data$characters <-
      c(
        list(
          model = model,
          shlib = libpath
        ),
        X$f$cgeneric$data$characters,
        Y$f$cgeneric$data$characters
      )

    ret$f$cgeneric$data$matrices <-
      c(
        X$f$cgeneric$data$matrices,
        Y$f$cgeneric$data$matrices
      )

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

    class(ret) <- "cgeneric"
    class(ret$f$cgeneric) <- "cgeneric"

    if(is.null(X$f$extraconstr)) {
      if(is.null(Y$f$extraconstr)) {
        if(debug)
          cat("No extraconstr!\n")
      } else {
        c2 <- Y$f$extraconstr
        ret$f$extraconstr <- list(
          A = kronecker(diag(X$f$n), c2$A),
          e = rep(c2$e, X$f$n)
        )
      }
    } else {
      c1 <- X$f$extraconstr
      if(is.null(Y$f$extraconstr)) {
        ret$f$extraconstr <- list(
          A = kronecker(c1$A, diag(Y$f$n)),
          e = rep(c1$e, each = Y$f$n)
        )
      } else {
        c2 <- Y$f$extraconstr
        ret$f$extraconstr <- list(
          A = rbind(
            kronecker(c1$A, diag(ncol(c2$A))),
            kronecker(diag(ncol(c1$A)), c2$A)
          ),
          e = c(rep(c1$e, each = ncol(c2$A)),
                rep(c2$e, ncol(c1$A)))
        )
      }
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
    rmodel <- list(
      f = list(
        model = "rgeneric",
        n = n,
        rgeneric = list(
          definition =
            compiler::cmpfun(
              kmodel,
              options = list(optimize = 3L)),
          debug = debug,
          optimize = TRUE
        )
      )
    )
    class(rmodel) <- "rgeneric"
    class(rmodel$f$rgeneric) <- "rgeneric"
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
    rmodel <- list(
      f = list(
        model = "rgeneric",
        n = n,
        rgeneric = list(
          definition =
            compiler::cmpfun(
              kmodel,
              options = list(optimize = 3L)),
          debug = debug,
          optimize = TRUE
        )
      )
    )
    class(rmodel) <- "rgeneric"
    class(rmodel$f$rgeneric) <- "rgeneric"
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

    return(INLA::inla.rgeneric.define(
      model = kmodel,
      optimize = TRUE
    ))

  }
)
