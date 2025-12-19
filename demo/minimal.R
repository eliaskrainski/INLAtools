suppressPackageStartupMessages({
	library(INLAtools)
	library(INLA)
})

R <- Sparse(crossprod(diff(diag(10))))
R

cmodel <- cgeneric("generic0", R = R, param = c(1, 0.01))

cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))

fit <- inla(formula = y ~ 0 + f(i, model = cmodel), 
	    data = data.frame(y = NA, i = 1:nrow(R)), 
	    control.family = cfam)

stopifnot(isTRUE(fit$ok))

