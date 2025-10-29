library(sf)

if(require(INLA)) {
    
library(INLAspacetime)

## define a graph on a grid
nbl <- spdep::grid2nb(d = c(7, 8))
graph.file <- "grid.graph"
spdep:::nb2INLA(file = graph.file, nb = nbl)

## build the spatial structure matrix
nn <- spdep::card(nbl)
n <- length(nn)
Gs <- sparseMatrix(
    i = rep(1:n, nn),
    j = unlist(nbl[nn>0]), x = 1L
)
Rs <- Diagonal(n = n, x = nn) - Gs

## simulate data
nt <- 50
stdf0 <- inla.knmodels.sample(
    graph.file, m = nt,
    tau.t = 2, tau.s = 2, tau.st = 3)

stdf <- data.frame(stdf0[c("time", "space", "spacetime")])
stdf$N <- (rgamma(n, 4, 1e-5)/10)[stdf$space]
stdf$x1 <- runif(n*nt)
stdf$x2 <- runif(n*nt)
stdf$x3 <- runif(n*nt)
beta.true <- c(-10, 0.5, -0.1, 0.03)
stdf$eta <- beta.true[1] + stdf0$x$eta +
    beta.true[2] * stdf$x1 +
    beta.true[3] * stdf$x2 +
    beta.true[4] * stdf$x3 
stdf$observed <- rpois(n*nt, stdf$N/(1+exp(-stdf$eta)))
stdf$expected <- stdf$N * (sum(stdf$observed)/sum(stdf$N))

c(sum(stdf$observed),  sum(stdf$expected),
  sum(stdf$observed) / sum(stdf$expected))

model00 <- observed ~ x1 + x2 + x3

## 'besagproper2' Kronecker with 'ar1'
model.p <- update(
    model00,
    .~.+ f(space, model = 'besagproper', 
           graph = graph.file, group = time,
           control.group = list(model = 'ar1')))
fit.p <- inla(
    formula = model.p,
    family = "poisson", E = expected,
    data = stdf,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

fit.p$cpu.used
fit.p$summary.fixed[, c(1,2,3,5)]
fit.p$summary.hyperpar[, c(1,2,3,5)]

## K-H type IV
tgraph <- sparseMatrix(i=c(2:nt, 1:(nt-1)), j=c(1:(nt-1), 2:nt), x=1)

model0 <- update(
    model00,
    .~.+f(space, model = 'bym2', graph = graph.file) +
        f(time, model = 'bym2', graph = tgraph))
fit.i4 <- inla.knmodels(
    formula = model0, data = stdf, 
    family = 'poisson', E = expected, 
    control.st = list(
        time = time, space = space, spacetime = spacetime,
        graph = graph.file, type = 4)
)

fit.i4$cpu.used
fit.i4$summary.fixed[, c(1,2,3,5)]
fit.i4$summary.hyper[, c(1,2,3,5)]


Rt <- Diagonal(n = nt, x = c(1, rep(2, nt-2), 1)) -tgraph
ct <- cgeneric("generic0", R = Rt, param = c(1, NA))

cs <- cgeneric("generic0", R = Rs, param = c(1, 0.05))

kh4 <- kronecker(ct, cs)

str(ct$f$extraconstr)
str(cs$f$extraconstr)
str(kh4$f$extraconstr)

nt+n
str(kh4$f$extraconstr)

image(prec(kh4, theta = 0))

fit4b <- inla(
    formula = update(model0, .~.+f(spacetime, model = kh4)), 
    data = stdf, family = 'poisson', E = expected
)

fit4b$cpu.used
fit4b$summary.fixed[, c(1,2,3,5)]
fit4b$summary.hyper[, c(1,2,3,5)]

grep("constraints", fit.i4$logfile, value = TRUE)
grep("constraints", fit4b$logfile, value = TRUE)

rbind(fit.i4$mode$theta, fit4b$mode$theta)



par(mfrow = c(3,4), mar = c(4, 4, 1, 1), mgp = c(3, 1, 0), las = 1, bty = "n")
for(k in 1:4) {
    plot(fit.p$marginals.fixed[[k]], type = 'l', lty = 2,
         xlab = names(fit.p$marginals.fixed)[k], ylab = "density")
    lines(fit.i4$marginals.fixed[[k]], lwd = 2)
    lines(fit4b$marginals.fixed[[k]], lty = 3, lwd = 2, col = 2)
}
for(k in 1:5) {
    plot(fit.i4$marginals.hyperpar[[k]], type = 'l', lwd = 2,
         xlab = names(fit.i4$marginals.hyperpar)[k], ylab = "density")
    if(k<5) {
        lines(fit4b$marginals.hyperpar[[k]], lty = 3, lwd = 2, col = 2)
    } else {
        lines(inla.tmarginal(exp, fit4b$marginals.hyperpar[[k]]),
              lty = 3, lwd = 2, col = 2)
    }
}
for(k in 1:3)
    plot(fit.p$marginals.hyperpar[[k]], type = 'l', lty = 2,
         xlab = names(fit.p$marginals.hyperpar)[k], ylab = "density")
    

rbind(fit.p$summary.hyperpar[1, , drop = FALSE],
      fit.i4$summary.hyperpar[5, , drop = FALSE],
      fit4b$summary.hyperpar[5, , drop = FALSE]
      )

cor(cbind(proper = fit.p$summary.random$space$mean,
          i4 = fit.i4$summary.random$spacetime$mean,
          k4 = fit.i4$summary.random$spacetime$mean))

}
