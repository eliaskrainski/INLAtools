library(INLAtools)

R <- Sparse(crossprod(diff(diag(10))))
R

m <- cgeneric("generic0", R = R,
	      scale = FALSE,
	      param = c(1, 0.01))
m

all.equal(R, Sparse(prec(m, theta = 0)))

graph(m)
prior(m, theta = 0)
prior(m, theta = matrix(-1:1, 1)) ## see ?prior.cgeneric

