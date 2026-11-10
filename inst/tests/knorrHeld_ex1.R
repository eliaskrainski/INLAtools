
setwd(here::here("vignettes"))

library(INLAtools)
library(sf)
library(spdep)
library(ggplot2)

data(GGHB.IZ, package = "CARBayesdata")
head(GGHB.IZ, 2)

ggplot(GGHB.IZ) + theme_minimal() + geom_sf()

data(pollutionhealthdata, package = "CARBayesdata")

(n <- nrow(GGHB.IZ))

head(pollutionhealthdata,2)

dataf <- data.frame(
    i = pmatch(pollutionhealthdata$IZ, GGHB.IZ$IZ, duplicates.ok = TRUE),
    t = pollutionhealthdata$year - min(pollutionhealthdata$year) + 1L,
    pollutionhealthdata[, c("pm10", "jsa", "price",
                            "observed", "expected")]    
)
head(dataf,2)

(m <- max(dataf$t))

## neighbour list
st_bbox(GGHB.IZ)

nbl0 <- poly2nb(GGHB.IZ)
nbl0
nbl <- poly2nb(GGHB.IZ, snap = 200)
nbl

par(mar = c(0,0,0,0))
plot(st_geometry(GGHB.IZ))
plot(nbl, st_coordinates(st_centroid(GGHB.IZ)), 
     add = TRUE, pch = 8, col = "red", lty = 2)
plot(nbl0, st_coordinates(st_centroid(GGHB.IZ)), 
     add = TRUE, pch = 0, col = "blue")
legend("bottom", c("0", "200m"), title = "snap", bty = "n",
       lty = 1:2, lwd = 1, col = c("blue", "red"))

(nn <- card(nbl))

G <- sparseMatrix(
    i = rep(1:n, nn),
    j = unlist(nbl[nn>0]),
    x = 1L,
    dims = c(n, n),
    repr = 'T'
)
R <- Diagonal(n = n, x = nn) - G

image(R)

###
H <- crossprod(diff(Diagonal(m), differences = 1))

if(n>10) {
    H[1:5, 1:5]
    H[-4:0+m, -4:0+m]
} else {
    H
}

## cgeneric models
cgt <- cgeneric(
    model = "generic0", R = R,
    constr = TRUE, scale = TRUE,
    param = c(1, 0.5)) ## P(sigma > 1) = 0.5
cgs <- cgeneric(
    model = "generic0", R = H,
    constr = TRUE, scale = TRUE,
    param = c(1, NA)) ## fix sigma, sigma = 1

## Kronecker of cgeneric models 1 and 2
cgst4 <- kronecker(cgt, cgs)

dim(cgt$f$extraconstr$A)
dim(cgs$f$extraconstr$A)
c(n,m)
dim(cgst4$f$extraconstr$A)

if(require(INLA)) {

## proper model
    m1f <- observed ~ jsa + price + pm10 + 
        f(i, model = "besagproper", graph = G, group = t,
          control.group = list(model = "ar1"))

    ccomp <- list(waic = TRUE, dic = TRUE, cpo = TRUE)
    fit1 <- inla(
        formula = m1f,
        data = dataf,
        family = "poisson", E = expected,
        control.compute = ccomp
    )

    round(fit1$summary.fixed[, c(1,2,3,5)], 4)
    round(fit1$summary.hyper[, c(1,2,3,5)], 4)

    ## (scaled) Knorr-Held type 4
    dataf$st <- 1:nrow(dataf)
    mst4 <- observed ~ jsa + price + pm10 +
        f(i, model = 'bym2', graph = G) +
        f(t, model = 'bym2', graph = H) +
        f(st, model = cgst4)

    fit4 <- inla(
        formula = mst4,
        data = dataf,
        family = "poisson", E = expected,
        control.compute = ccomp
    )

    1 + 1 + (n*m)-(n-1)*(m-1)
    grep("onstraints", fit4$logfile, value = TRUE)

    gstats <- function(x) {
        c(waic=x$waic$waic, dic = x$dic$dic,
          cpo = -sum(log(x$cpo$cpo)))
    }

    gcv1 <- inla.group.cv(fit1, num.level.sets = 3)
    gcv4 <- inla.group.cv(fit4, num.level.sets = 3)

    summary(sapply(gcv1$groups, function(x) length(x$idx)))
    summary(sapply(gcv4$groups, function(x) length(x$idx)))

    summary(sapply(gcv1$groups, function(x) mean(x$corr[-1])))
    summary(sapply(gcv4$groups, function(x) mean(x$corr[-1])))

    rbind(m1 = c(gstats(fit1), gcv=-sum(log(gcv1$cv))),
          m4 = c(gstats(fit4), gcv=-sum(log(gcv4$cv))))

    (sfix <- round(fit4$summary.fixed[, c(1,2,3,5)], 4))
    (fnams <- c("0", rownames(sfix)[-1]))
    
    round(fit4$summary.hyperpar[, c(1,2,3,5)], 4)
    
    par(mfrow = c(3,3), mar = c(3,3,0.5,0.5), mgp = c(2,0.5,0), bty = 'n')
    for(k in 1:4)
        plot(inla.smarginal(fit4$marginals.fixed[[k]]), type = "l",
             xlab = as.expression(bquote(beta[.(fnams[k])])))
    plot(inla.tmarginal(function(x) exp(-x/2),
                        fit4$internal.marginals.hyperpar[[1]]), type = "l",
         xlab = expression(sigma[t]))
    plot(inla.smarginal(fit4$marginals.hyperpar[[2]]), type = "l",
         xlab = expression(psi[t]))
    plot(inla.tmarginal(function(x) exp(-x/2),
                        fit4$internal.marginals.hyperpar[[3]]), type = "l",
         xlab = expression(sigma[s]))
    plot(inla.smarginal(fit4$marginals.hyperpar[[4]]), type = "l",
         xlab = expression(psi[t]))
    plot(inla.tmarginal(function(x) exp(x),
                        fit4$internal.marginals.hyperpar[[5]]), type = "l",
         xlab = expression(sigma[s]))
    
}

