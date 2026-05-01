library(INLAtools)

R <- Sparse(crossprod(diff(diag(10))))
R

m <- cgeneric("generic0", R = R,
	      scale = FALSE,
	      param = c(1, 0.01))
m

all.equal(R, Sparse(cgeneric_Q(m, theta = 0)))

cgeneric_graph(m)
cgeneric_prior(m, theta = 0)
cgeneric_prior(m, theta = matrix(-1:1, 1)) ## see ?prior.cgeneric

