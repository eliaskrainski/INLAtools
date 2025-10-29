library(INLAtools)

if(require(INLA)) {
    
### define a graph on a grid
    nbl <- spdep::grid2nb(d = c(5, 10))
    graph.file <- "grid.graph"
    spdep:::nb2INLA(file = graph.file, nb = nbl)
    
### build the spatial structure matrix
    nn <- spdep::card(nbl)
    n <- length(nn)
    Gs <- sparseMatrix(
        i = rep(1:n, nn),
        j = unlist(nbl[nn>0]), x = 1L
    )
    Rs <- Diagonal(n = n, x = nn) - Gs
    
### simulate latent fields from the (scaled) Knorr-Held type IV
    nt <- 30
    stdf0 <- inla.knmodels.sample(
        graph.file, m = nt,
        tau.t = 2, tau.s = 2, tau.st = 3)

### simulate data (counts)
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

### 'besagproper2' Kronecker with 'ar1'
    model.prop <- update(
        model00,
        .~. + f(space, model = 'besagproper', 
                graph = graph.file, group = time,
                control.group = list(model = 'ar1')))
    fit.prop <- inla(
        formula = model.prop,
        family = "poisson", E = expected,
        data = stdf,
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
    
    fit.prop$cpu.used
    fit.prop$summary.fixed[, c(1,2,3,5)]
    fit.prop$summary.hyperpar[, c(1,2,3,5)]
    
### K-H type IV
    tgraph <- sparseMatrix(i=c(2:nt, 1:(nt-1)),
                           j=c(1:(nt-1), 2:nt), x=1)
    model0 <- update(
        model00,
        .~.+f(space, model = 'bym2', graph = graph.file) +
            f(time, model = 'bym2', graph = tgraph))
    fit.kh4 <- inla.knmodels(
        formula = model0, data = stdf, 
        family = 'poisson', E = expected, 
        control.st = list(
            time = time, space = space, spacetime = spacetime,
            graph = graph.file, type = 4),
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
    
    fit.kh4$cpu.used
    fit.kh4$summary.fixed[, c(1,2,3,5)]
    fit.kh4$summary.hyper[, c(1,2,3,5)]
    
    
    Rt <- Diagonal(n = nt, x = c(1, rep(2, nt-2), 1)) -tgraph
    ct <- cgeneric("generic0", R = Rt, param = c(1, NA))
    
    cs <- cgeneric("generic0", R = Rs, param = c(1, 0.05))
    
    kh4cg <- kronecker(ct, cs)
    
    str(ct$f$extraconstr)
    str(cs$f$extraconstr)
    str(kh4cg$f$extraconstr)
    
    nt+n
    str(kh4cg$f$extraconstr)
    
    image(prec(kh4cg, theta = 0))
    
    fit.kh4cg <- inla(
        formula = update(model0, .~.+f(spacetime, model = kh4cg)), 
        data = stdf, family = 'poisson', E = expected,
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
    
    fit.kh4cg$cpu.used
    fit.kh4cg$summary.fixed[, c(1,2,3,5)]
    fit.kh4cg$summary.hyper[, c(1,2,3,5)]
    
    grep("constraints", fit.kh4$logfile, value = TRUE)
    grep("constraints", fit.kh4cg$logfile, value = TRUE)

    grep("fn-calls", fit.kh4$logfile, value = TRUE)
    grep("fn-calls", fit.kh4cg$logfile, value = TRUE)
    
    rbind(fit.kh4$mode$theta, fit.kh4cg$mode$theta)
    
    par(mfrow = c(3,4), mar = c(4, 4, 1, 1), mgp = c(3, 1, 0), las = 1, bty = "n")
    for(k in 1:4) {
        plot(fit.prop$marginals.fixed[[k]], type = 'l', lty = 2,
             xlab = names(fit.prop$marginals.fixed)[k], ylab = "density")
        lines(fit.kh4$marginals.fixed[[k]], lwd = 2)
        lines(fit.kh4cg$marginals.fixed[[k]], lty = 3, lwd = 2, col = 2)
    }
    for(k in 1:5) {
        plot(fit.kh4$marginals.hyperpar[[k]], type = 'l', lwd = 2,
             xlab = names(fit.kh4$marginals.hyperpar)[k], ylab = "density")
        if(k<5) {
            lines(fit.kh4cg$marginals.hyperpar[[k]], lty = 3, lwd = 2, col = 2)
        } else {
            lines(inla.tmarginal(exp, fit.kh4cg$marginals.hyperpar[[k]]),
                  lty = 3, lwd = 2, col = 2)
        }
    }
    for(k in 1:3)
        plot(fit.prop$marginals.hyperpar[[k]], type = 'l', lty = 2,
             xlab = names(fit.prop$marginals.hyperpar)[k], ylab = "density")
    
    
    rbind(fit.prop$summary.hyperpar[1, , drop = FALSE],
          fit.kh4$summary.hyperpar[5, , drop = FALSE],
          fit.kh4cg$summary.hyperpar[5, , drop = FALSE]
          )
    
    cor(cbind(proper = fit.prop$summary.random$space$mean,
              kh4 = fit.kh4$summary.random$spacetime$mean,
              k4 = fit.kh4$summary.random$spacetime$mean))
    
    gcv <- lapply(list(proper = fit.prop, kh4 = fit.kh4, kh4cg = fit.kh4cg),
                  inla.group.cv, num.level.sets = 10)

### compare some statistics
    rbind(mlik = c(proper = fit.prop$mlik[1],
                   kh4 = fit.kh4$mlik[1],
                   kh4 = fit.kh4cg$mlik[1]),
          dic = c(fit.prop$dic$dic,
                  fit.kh4$dic$dic,
                  fit.kh4cg$dic$dic),
          waic = c(fit.prop$waic$waic,
                   fit.kh4$waic$waic,
                   fit.kh4cg$waic$waic),
          cpo = c(-sum(log(fit.prop$cpo$cpo)),
                  -sum(log(fit.kh4$cpo$cpo)),
                  -sum(log(fit.kh4cg$cpo$cpo))),
          gcv = sapply(gcv, function(x) -sum(log(x$cv))))

}
