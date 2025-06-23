library(INLAtools)

## first dim
(n1 <- nrow(
     G1 <- sparseMatrix(
         i = c(2, 3, 1, 4, 1, 4, 5, 2, 3, 3),
         j = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5),
         )
 ))
R1 <- Diagonal(n = n1, x = colSums(G1)) - G1
R1

## 2nd dim
(n2 <- nrow(
     R2 <- sparseMatrix(
         i = c(1L, 1L, 2L, 2L, 2L, 3L, 3L),
         j = c(1L, 2L, 1L, 2L, 3L, 3L, 2L),
         x = c(2,-1, -1,3,-1, 4, -1))
 ))
R2

R12 <- kronecker(R1, R2)
R12

## cgeneric models
cg1 <- cgeneric(
    model = "generic0", R = R1,
    constr = FALSE, scale = FALSE,
    param = c(1, 0.5)) ## P(sigma > 1) = 0.5
cg2 <- cgeneric(
    model = "generic0", R = R2,
    constr = FALSE, scale = FALSE,
    param = c(1, NA)) ## fix sigma, simga = 1

## Kronecker of cgeneric models 1 and 2
cg12 <- kronecker(cg1, cg2)

all.equal(Sparse(R12),
          Sparse(prec(cg12, theta = 0)))

if("INLA" %in% loadedNamespaces()) {

    ## create fake data to call inla()
    data2 <- as.data.frame(
        expand.grid(i1 = 1:n1,
                    i2 = 1:n2)
    )
    ## additional index (combined i1 and i2)
    data2$ii <- 1:nrow(data2)
    ## no data
    data2$y <- rep(NA, nrow(data2))

    (n1*n2)==nrow(data2)

    R1
    head(data2,5)

    m1f <- y ~ 0 +
        f(ii, model = 'generic0', Cmatrix = R12)

    hfix <- list(prec = list(initial = 10, fixed = TRUE))

    fit1 <- INLA::inla(
        formula = m1f,
        data = data2,
        control.mode = list(theta = 0, fixed = TRUE),
        control.family = list(hyper = hfix),
        control.compute = list(config = TRUE)
    )

    all.equal(Sparse(R12),
              Sparse(prec(fit1)))

    m2f <- y ~ 0 +
        f(i2, model = 'generic0', Cmatrix = R2,
          group = i1,
          control.group = list(
              model = 'besag', graph = G1,
              scale.model = FALSE))

    fit2 <- INLA::inla(
        formula = m2f,
        data = data2,
        control.mode = list(theta = 0, fixed = TRUE),
        control.family = list(hyper = hfix),
        control.compute = list(config = TRUE)
    )

    all.equal(Sparse(R12),
              Sparse(prec(fit2)))

}
