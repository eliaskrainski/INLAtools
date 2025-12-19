library(INLAtools)

R <- Sparse(crossprod(diff(diag(10))))
R

m <- cgeneric("generic0", R = R, scale = FALSE, param = c(1, 0.01))

all.equal(R, Sparse(prec(m, theta = 0)))

## see ?prior for more

