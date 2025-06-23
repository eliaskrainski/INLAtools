
library(INLAtools)

## first dim
G1 <- sparseMatrix(
    i = c(2, 3, 1, 4, 1, 4, 2, 3),
    j = c(1, 1, 2, 2, 3, 3, 4, 4),
)
R1 <- Diagonal(n = nrow(G1), x = colSums(G1)) - 0.9*G1
R1

## 2nd dim
R2 <- sparseMatrix(
    i = c(1L, 1L, 2L, 2L, 2L, 3L, 3L),
    j = c(1L, 2L, 1L, 2L, 3L, 3L, 2L),
    x = c(2,-1, -1,3,-1, 4, -1))
R2

(r <- 0.9)
R3 <- sparseMatrix(
    i = c(1:3, 1:2, 2:3),
    j = c(1:3, 2:3, 1:2),
    x = c(1, 1+r^2, 1, rep(-r, 4))
)/(1-r^2)
R3

c(n1 <- nrow(R1),
  n2 <- nrow(R2),
  n3 <- nrow(R3))

### cgeneric models for each dim
cg1 <- cgeneric(
    model = "generic0", R = R1,
    constr = FALSE, scale = FALSE,
    param = c(1, 0.5))
cg2 <- cgeneric(
    model = "generic0", R = R2,
    constr = FALSE, scale = FALSE,
    param = c(1, NA)) ## fixed to 1
cg3 <- cgeneric(
    model = "generic0", R = R3,
    constr = FALSE, scale = FALSE,
    param = c(1, NA)) ## fixed to 1

## kronecker 2x1
cg21 <- kronecker(cg2, cg1)
cg21Q <- prec(cg21, theta = 0)

## compare
R21 <- kronecker(R2, R1)
all.equal(Sparse(R21),
          Sparse(cg21Q))

## kronecker 3x(2x1)
cg321 <- kronecker(cg3, cg21)
cg321Q <- prec(cg321, theta = 0)

R321 <- kronecker(R3, R21)
all.equal(Sparse(R321),
          Sparse(cg321Q))

if("INLA" %in% loadedNamespaces()) {

    dataf <- as.data.frame(
    expand.grid(i1 = 1:n1,
                i2 = 1:n2,
                i3 = 1:n3)
    )
    nrow(dataf)

    ## additional indexes (combined i1 and i2), and 1:(n1*n2*n3)
    dataf$i12 <- rep(1:(n1*n2), n3)
    dataf$iii <- 1:nrow(dataf)

    ## no data
    dataf$y <- rep(NA, nrow(dataf))

    ## use the group
    m1f <- y ~ 0 +
        f(i12, model = 'generic0', Cmatrix = R21,
          group = i3,
          control.group = list(
              model = 'ar1',
              hyper = list(rho = list(
                               initial = log((1+r)/(1-r)),
                               fixed = TRUE))))

    hfix <- list(prec = list(initial = 10, fixed = TRUE))

    fit1 <- INLA::inla(
        formula = m1f,
        data = dataf,
        control.mode = list(theta = 0, fixed = TRUE),
        control.family = list(hyper = hfix),
        control.compute = list(config = TRUE)
    )

    Q1 <- prec(fit1)

    all.equal(Sparse(R321),
              Sparse(Q1))

    ## 'fit' cg321
    m321 <- y ~ 0 + f(iii, model=cg321)
    fit2 <- INLA::inla(
        formula = m321,
        data = dataf,
        control.mode = list(theta = 0, fixed = TRUE),
        control.family = list(hyper = hfix),
        control.compute = list(config = TRUE)
    )


    Q2 <- prec(fit2)

    all.equal(Sparse(R321),
              Sparse(Q2))

}
