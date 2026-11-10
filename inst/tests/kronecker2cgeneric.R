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

R12 <- kronecker(Sparse(R1), Sparse(R2))
R12

## build a cgeneric generic0 model with R1, Q = \tau * R1: model 1
cg1 <- cgeneric(
    model = "generic0", R = R1,
    constr = FALSE, scale = FALSE,
    param = c(1, 0.5)) ## P(sigma > 1) = 0.5


## build a cgeneric generic0 model with R2,
## but fix the precision to 1, Q = 1 * R2: model 2 
cg2 <- cgeneric(
    model = "generic0", R = R2,
    constr = FALSE, scale = FALSE,
    param = c(1, NA)) ## fix sigma, simga = 1

## build a cgeneric model where the precision is 
## the Kronecker Q = \tau R1 (o) R2 
cg12 <- kronecker(cg1, cg2)

## Check
stopifnot(all.equal(
    Sparse(R12),
    Sparse(cgeneric_Q(cg12, theta = 0))))

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
    
    Qfitted <- fit0$misc$configs$config[[1]]$Q
    stopifnot(all.equal(Sparse(upperPadding(R12)),
                        Sparse(Qfitted)))
    
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

    Qfitted1 <- fit1$misc$configs$config[[1]]$Q
    stopifnot(all.equal(Sparse(upperPadding(R12)),
                        Sparse(Qfitted1)))

}

