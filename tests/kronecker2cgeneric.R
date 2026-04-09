### Q = \tau (R2 \otimes R1)

library(INLAtools)

## first
(n1 <- nrow(
     G1 <- sparseMatrix(
         i = c(2, 3, 1, 4, 1, 5, 2, 3),
         j = c(1, 1, 2, 2, 3, 3, 4, 5),
         )
 ))
R1 <- Diagonal(n = n1, x = colSums(G1)) - G1
R1

## 2nd 
(n2 <- nrow(
     G2 <- sparseMatrix(
         i = c(1L, 2L, 2L, 3L),
         j = c(2L, 1L, 3L, 2L),
         dims = c(4L, 4L)
     )
 )
)
R2 <- Diagonal(n = n2, x = colSums(G2)) - G2
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

## Checks
g1 <- graph(cg1, optimize = !TRUE)
g2 <- graph(cg2, optimize = !TRUE)

g1
g2

if(require(INLA)) {

    ## create fake data to call inla()
    data2 <- data.frame(
        expand.grid(i1 = 1:n1,
                    i2 = 1:n2),
        y = NA
    )

    mfg <- y ~ 0 +
        f(i2, model = 'generic0', Cmatrix = R2,
          group = i1,
          control.group = list(
              model = 'besag', graph = G1,
              scale.model = FALSE))

    hfix <- list(prec = list(initial = 10, fixed = TRUE))

    fit0 <- inla(
        formula = mfg,
        data = data2,
        control.mode = list(theta = 0, fixed = TRUE),
        control.family = list(hyper = hfix),
        control.compute = list(config = TRUE)
    )
    
    print(all.equal(Sparse(R12),
                    Sparse(prec(fit0))))

    ## overall index 
    (n1*n2)==nrow(data2)
    data2$ii <- 1:nrow(data2)

    mfcgk <- y ~ 0 +
        f(ii, model = 'generic0', Cmatrix = R12)

    fit1 <- inla(
        formula = mfcgk,
        data = data2,
        control.mode = list(theta = 0, fixed = TRUE),
        control.family = list(hyper = hfix),
        control.compute = list(config = TRUE)
    )

    print(
        all.equal(Sparse(R12),
                  Sparse(prec(fit1)))
    )

}

