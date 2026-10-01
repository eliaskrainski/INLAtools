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
#' @param id.param integer to assign with param of the Mc model
#' would be ordered. Used if Mc has more than one parameter.
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
           id.param,
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
      minimum_version = "26.10.1")
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

    ## initial to count the number of parameters
    iniMc <- cgeneric_get(
      model = model,
      cmd = "initial")
    nthMc <- length(iniMc)
    if((missing(id.param)) || is.null(id.param)) {
      id.param <- min(1:nthMc)
    }
    stopifnot(id.param %in% 1:nthMc) ## work if nthMc==0

    ## Q1 graph
    graphMc <- cgeneric_get(
      model = model,
      cmd = "graph",
      optimize = FALSE)
    ## order it (row wise) as needed later
    oiMc <- order(graphMc@i)
    graphMc@i <- graphMc@i[oiMc]
    graphMc@j <- graphMc@j[oiMc]
    ugMc <- upperPadding(graphMc)
    ## indexing
    McU <- length(ugMc@i)
    idxLowfn <- function(i,j,n,m) {
      il <- which(i<j)
      nl <- length(il)
      idx <- 1:m
      if(nl>0) {
        r <- c(n*i+j, j[il]*n+i[il])
        o <- order(r)
        idx <- c(idx, idx[il])[o]
      }
      return(idx)
    }
    idxMc <- idxLowfn(ugMc@i, ugMc@j, nMc, McU)

    McAll <- McU*2-nMc
    stopifnot(McAll==length(idxMc))

    uGraph <- upperPadding(
      kronecker(matrix(1,K,K), graphMc))
    M <- McU*K + (McAll * (K*(K-1)/2))
    stopifnot(M==length(uGraph@i))
    if(dotArgs$debug) {
      cat('McU:', McU, "McAll:", McAll, "M:", M, "\n")
    }

    ## output order
    ## internally it compute Q by blocks of size nxn by row
    ## the diagonal blocks have McU elements
    ## the upper diagonal blocks have McAll elements
    k <- 0
    id0 <- vector("list", K)
    for(i in 1:(K-1)) {
      ok <- order(c(ugMc@i, rep(graphMc@i, K-i)))
      id0[[i]] <- k + ok
      k <- max(id0[[i]])
    }
    id0[[K]] <- k + 1:McU

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
            ## inla.cgeneric is needed to support INLA before August 2025
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
      c(list(n=as.integer(
        c(nMc=nMc, ## 0
          ndataMc,  ## 1:5
          Mc=McU, ## 6
          K, ## 7
          c(i=3,d=3,c=2,m=0,s=0), ## 8:12
          N=N, ## 13
          M=M ## 14
        ))),
        model$f$cgeneric$data$ints[-1],
        list(idParam = as.integer(id.param-1L),
             idxMc = as.integer(idxMc-1L),
             ii = uGraph@i,
             jj = uGraph@j,
             order = as.integer(unlist(id0)-1)
             )
        )

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
      c(list(model = cmodel,
             shlib = dotArgs$shlib),
        model$f$cgeneric$data$characters)

    ## matrices
    if(ndataMc[4]>0) {
      retModel$f$cgeneric$data$matrices <-
        model$f$cgeneric$data$matrices
    }

    ## sparse matrices (none added)
    retModel$f$cgeneric$data$smatrices <- model$f$cgeneric$data$smatrices

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
