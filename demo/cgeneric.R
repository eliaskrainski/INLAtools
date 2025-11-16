R <- Sparse(crossprod(diff(diag(10))))
m <- cgeneric("generic0", R = R, param = c(1, 0.01))
all.equal(R, prec(m))

## see ?prior for more

