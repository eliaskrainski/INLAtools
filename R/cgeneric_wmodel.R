#' Build a `cgeneric` object for a `wmodel`.
#' @description
#' Build data to implement the K dimensional W-model
#' applied to a given model, see details.
#' @param model a `cgeneric` model, the model 'Mc'.
#' @param K dimension of W.
#' @param lambda the parameter for the exponential prior on
#' the radius of the sphere, see details.
#' This can be a constant or vector up to length K.
#' If a length K `lambda` is provided, each one will be used to
#' define the PC-prior for each row-vector in W.
#' @param sigma.prior.reference numeric vector to set the reference
#' for each overall standard deviation parameter for its PC-prior.
#' Use this if the given `cgeneric` model 'Mc' has fixed variance.
#' @param sigma.prior.probability numeric vector with to
#' set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
#' This is the default case, that can be used when assuming 'Mc'
#' with unknown variance.
#' @param ... possible additional arguments,
#' such as debug, useINLAprecomp, shlib,
#' passed on to `cgeneric`.
#' @details
#' The precision matrix for the 'Mc' model is defined as
#'  \deqn{Q = (W (o) I_n) bdiag(Q_1, ..., Q_K) (W (o) I_n)'}
#' where the matrices \eqn{Q_j}, \eqn{j=1,...,K} are build from
#' a given model. The special case when \eqn{Q_j} is from a
#' CAR model is detailed in Martinez-Beneito, M.A., (2013).
#' The W matrix is a K by K matrix whose rows domain is the
#' surface of a K-1 sphere with unit radius, by setting
#' \eqn{\sum_j W_{ij}^2 = 1}, which imposes a constraint so that
#' there is K-1 unknown in each row-vector in W.
#' Additionally, we impose \eqn{W_{ii}>0}.
#' W is parametrized from \eqn{m=K(K-1} internal parameters.
#'  \eqn{\theata_0, \ldots, \theta_m},
#' where the first K=1 parameters are to define  the first row of W,
#' the next K-1 to the second row and so on. The map from the
#' internal parameters is done from the K-1 internal parameters to
#' each row of W by the three steps:
#' 1) set \eqn{x_{ii}=1}; 2) set \eqn{x_{-ii}=\theta}; and
#' 3) compute \eqn{W_{ij}=x_{ij}/\sqrt{\sum_j x_{ij}^2}}.
#' @references
#' Martinez-Beneito, M.A., (2013).
#' A general modelling framework for multivariate disease mapping.
#' Biometrika 100, 539–553. <doi:10.1093/biomet/ast023>.
#' @return a `cgeneric` object, see [cgeneric-class()].
#' @export
cgeneric_wmodel <-
  function(model,
           K,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           ...) {

    dotArgs <- list(...)
    if(is.null(dotArgs$debug)) {
      dotArgs$debug <- FALSE
    }
    K <- as.integer(K)[1]
    stopifnot(K>2)
    K2 <- K * K

    stopifnot(inherits(model, "cgeneric"))
    if(dotArgs$debug) {
      print(model)
    }
    stopifnot((nMc <- as.integer(model$f$n))>1)
    N <- K * nMc

    if(dotArgs$debug) {
      print(c(K=K, nMc=nMc, N = N))
    }

    stopifnot(all(lambda>0))
    if(length(lambda)==1) {
      lambda <- rep(lambda, K)
    }
    if(length(lambda)!=K) {
      stop('"length(lambda)" shold be 1 or K!')
    }

    pcSigmas <- pcParamCheck(
      npars = K,
      reference = sigma.prior.reference,
      probability = sigma.prior.probability
    )
    if(dotArgs$debug) {
      print(str(list(pcSigmas=pcSigmas)))
    }

    if(is.null(dotArgs$useINLAprecomp)) {
      dotArgs$useINLAprecomp <- TRUE
    }
    INLAvcheck <- packageCheck(
      name = "INLA",
      minimum_version = "26.9.24")
    if(is.na(INLAvcheck) & dotArgs$useINLAprecomp) {
      dotArgs$useINLAprecomp <- FALSE
      warning("INLA version is old. Setting 'useINLAprecomp = FALSE'!")
      warning("Using the developing internal model!")
      dotArgs$developing <- TRUE
    }
    cmodel <- "inla_cgeneric_wmodel"
    if(!is.null(dotArgs$developing)) {
      cmodel <- paste0(cmodel, "_dev")
    }

    dotArgs$shlib <- cgeneric_shlib_path(
      package = "INLAtools",
      useINLAprecomp = dotArgs$useINLAprecomp,
      debug = dotArgs$debug)

    ## Q1 graph and indexing
    graphMc <- cgeneric_get(
      model = model,
      cmd = "graph",
      optimize = FALSE)
    w1I <- kronecker(matrix(1,K,K), Diagonal(n=nMc))
    Gkk <- bdiag(lapply(1:K, function(k) graphMc))
    Graph <- Sparse(w1I %*% Gkk %*% Matrix::t(w1I))

    ugMc <- upperPadding(graphMc)
    McU <- length(ugMc@x)
    ugMc@x <- as.numeric(1:McU)
    graphMc@x <- Sparse(ugMc + t(ugMc) -
                          Diagonal(nMc, diag(ugMc)))@x

    Mc <- length(graphMc@x)
    M <- McU*K + (Mc * (K*(K-1)/2))
    if(dotArgs$debug) {
      cat('McU:', McU, "Mc:", Mc, "M:", M, "\n")
      print(list(graphMc=graphMc,
                 ugraphMc=ugMc,
                 Graph = Graph))
    }

    ## initial data
    retModel <- structure(
      list(
        f = list(
          model = "cgeneric",
          n = as.integer(N),
          cgeneric = structure(
            list(
              model = cmodel,
              shlib = dotArgs$shlib,
              n = as.integer(N),
              debug = as.integer(dotArgs$debug)
            ),
            # inla.cgeneric is needed to support INLA before August 2025
            class = c("inla.cgeneric.f", "inla.cgeneric")
          )
        )
      ),
      class = c("cgeneric", "inla.cgeneric")
    )
    retModel$f$cgeneric$data <- vector("list", 5L)
    names(retModel$f$cgeneric$data) <- c(
      "ints", "doubles", "characters",
      "matrices", "smatrices")

    ndataMc <- sapply(model$f$cgeneric$data, length)

    ## ints
    retModel$f$cgeneric$data$ints <-
      c(list(n=ndataMc), model$f$cgeneric$data$ints[-1], K=K)

    ## doubles
    if(ndataMc[2]>0) {
      retModel$f$cgeneric$data$doubles <-
        c(model$f$cgeneric$data$doubles,
          list(lambda = as.double(lambda)))
    } else {
      retModel$f$cgeneric$data$doubles <-
        list(lambda = as.double(lambda))
    }

    ## characters
    retModel$f$cgeneric$data$characters <-
      model$f$cgeneric$data$characters

    ## matrices
    if(ndataMc[4]>0) {
      retModel$f$cgeneric$data$matrices <-
        model$f$cgeneric$data$matrices
    }

    if(ndataMc[5]>0) {
      retModel$f$cgeneric$data$smatrices <-
        c(model$f$cgeneric$data$smatrices, list(
          Qgraph = c(nrow(Graph), ncol(Graph), length(Graph@x),
                     Graph@i, Graph@j, Graph@x)))
    } else {
      retModel$f$cgeneric$data$smatrices <-
        list(Qgraph = c(nrow(Graph), ncol(Graph), length(Graph@x),
                        Graph@i, Graph@j, Graph@x))
    }

## setup extraconstr
    Ae <- dotArgs$extraconstr
    if(!is.null(Ae)) {
      stopifnot(all(c("A", "e") %in% names(Ae)))
      stopifnot(ncol(Ae$A)==nMc)
      stopifnot(nrow(Ae$A)==length(Ae$e))
    }
    if((!is.null(Ae)) | (!is.null(model$f$extraconstr))) {
       retModel$f$extraconstr <-
        kronecker_extraconstr(
          Ae, model$f$extraconstr, K, nMc)
      if(dotArgs$debug) {
        cat(nrow(rmodel$f$extraconstr),
            " 'extraconstr' built\n")
      }
    }

    ## this is for inlabru mapper system to work
    ## Note: kron(X,Y) gives X-major ordering (the X-index varies slowly, the
    ## Y-index varies quickly), which requires the mappers to be in reverse
    ## order, multi(Y,X):
    retModel$mapper <- multi_generic_model_mapper(
      list(model, mapper1(list(n=K))))

    return(retModel)

  }
