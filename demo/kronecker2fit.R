library(INLAtools)

## first model precision structure
n1 <- 15
(R1 <- crossprod(diff(Diagonal(n1))))


## second model precision structure 
G2 <- sparseMatrix(
  i = c(1,1,2,2,3,3,6), 
  j = c(2,3,4,6,4,5,7), 
  symmetric = TRUE)
(n2 <- nrow(G2)) 
R2 <- Diagonal(n = n2, x = colSums(G2)) - G2
R2

## make cgeneric for each model
## cgeneric model 1
cg1 <- cgeneric(
    model = "generic0", 
    R = R1, ## precision structure matrix
    param = c(1, 0.05)) ## P(sigma > 1) = 0.05


## cgeneric model 2
cg2 <- cgeneric(
    model = "generic0", 
    R = R2, ## precision structure matrix
    param = c(1, NA)) ## fix sigma = 1


## cgeneric Kronecker model
cg12 <- kronecker(cg1, cg2)

## inspect
(N <- cg12$f$n)
tau <- 4
Q <- cgeneric_Q(cg12, theta = log(tau))

image(Q)


## model fitting
if(require(INLA)) {
    
    k <- 5
    xx <- inla.qsample(
        n = k, 
        Q = Q + Diagonal(N, 1e-9), 
        constr = cg12$f$extraconstr, 
        seed = 1
    )
    apply(xx, 2, summary)
    
    ## organize the data
    dataf <- data.frame(
        i1 = rep(rep(1:n1, each = n2), k),
        i2 = factor(rep(rep(1:n2, n1), k)),
        i = rep(1:N, k), 
        r = rep(1:k, each = N),
        x = as.vector(xx)
    )
    head(dataf, 10)

    ## ggplot 
    library(ggplot2)
    ggplot(dataf) + theme_minimal() + 
        geom_line(aes(x = i1, y = x, group = i2, color = i2)) + 
        facet_wrap(~r)
    
    ## the outcome 
    dataf$y <- rpois(N * k, exp(3 + dataf$x))

    ## inla call
    fit <- inla(
        formula = y ~ f(i, model = cg12, replicate = r),
        family = "poisson", 
        data = dataf)
    
    ## summary of posterior marginals 
    fit$summary.fixed
    fit$summary.hyperpar

    ## randon effects posterior mean aganist simulated
    par(mar = c(3,3,0.5,0.5), mgp = c(2,0.5,0), bty = "n")
    plot(fit$summary.random$i$mean, xx, pch = 8,
         xlab = "Posterior mean", ylab = "Simulated")

    ## sigma posterior marginal 
    pm.sigma <- inla.tmarginal(
        fun = function(x) exp(-x/2), ## from log(tau) to sigma
        marginal = fit$internal.marginals.hyperpar[[1]])
    
    1/sqrt(tau)
    inla.zmarginal(pm.sigma)

    plot(pm.sigma, type = 'l',
         xlab = expression(sigma),
         ylab = as.expression(bquote(P(sigma~"|"~y))))

}
