library(INLAtools)

### Q = \tau (R2 \otimes R1)

## first
(n1 <- nrow(
     G1 <- sparseMatrix(
         i = c(2, 3, 1, 4, 1, 4, 5, 2, 3, 3),
         j = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5),
         )
 ))
R1 <- Diagonal(n = n1, x = colSums(G1)) - G1
R1

## 2nd 
(n2 <- nrow(
     R2 <- sparseMatrix(
         i = c(1L, 1L, 2L, 2L, 2L, 3L, 3L),
         j = c(1L, 2L, 1L, 2L, 3L, 3L, 2L),
         x = c(2,-1, -1,3,-1, 4, -1))
 ))
R2

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

all.equal(Sparse(kronecker(R1, R2)),
          Sparse(prec(cg12, theta = 0)))



