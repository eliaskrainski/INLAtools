## Q = \tau (Rt \otimes Rs)
##  Rt is the temporal structure matrix
##  Rs is the spatial structure matrix

library(sf)
library(spdep)
library(Matrix)
library(rnaturalearth)
library(INLAtools)

map <- ne_states(country = "Italy") ## "United States of America")
nrow(map)

plot(st_geometry(map))

nbl <- poly2nb(map)
nn <- card(nbl)
n <- length(nn)

G <- sparseMatrix(
    i = rep(1:n, nn),
    j = unlist(nbl)[nn>0],
    x = 1L
)
Rs <- Diagonal(n = n, x = nn) - G

###
m <- 30
Rt <- crossprod(diff(Diagonal(m)))
Rt[1:3, 1:3]
Rt[-2:0+m, -2:0+m]

Rst <- kronecker(Rt, Rs)
image(Qst)

## cgeneric models
cgt <- cgeneric(
    model = "generic0", R = Rt,
    constr = FALSE, scale = FALSE,
    param = c(1, 0.5)) ## P(sigma > 1) = 0.5
cgs <- cgeneric(
    model = "generic0", R = Rs,
    constr = FALSE, scale = FALSE,
    param = c(1, NA)) ## fix sigma, sigma = 1

## Kronecker of cgeneric models 1 and 2
cgst <- kronecker(cgt, cgs)

all.equal(Sparse(Rst),
          Sparse(prec(cgst, theta = 0)))

if(require(INLA)) {

    ## create fake data to call inla()
    dataf <- data.frame(
        t = rep(1:m, each = n),
        s = rep(1:n, m),
        i = 1:(n*m),
        y = NA
    )
    
    m1f <- y ~ 0 +
        f(i, model = 'generic0', Cmatrix = Rst)

    fit1 <- inla(
        formula = m1f,
        data = dataf,
        control.mode = list(theta = c(0,0), fixed = TRUE)
    )

    all.equal(Sparse(Rst),
              Sparse(prec(fit1)))

    m2f <- y ~ 0 +
        f(s, model = 'generic0', Cmatrix = Rs,
          group = t,
          control.group = list(
              model = 'rw1',
              scale.model = FALSE))

    fit2 <- inla(
        formula = m2f,
        data = dataf,
        control.mode = list(theta = c(0,0), fixed = TRUE)
    )

    print(
        all.equal(Sparse(Rst),
                  Sparse(prec(fit2)))
    )

    fit3 <- inla(
        y ~ 0 + f(i, model = cgst), data = dataf,
        control.mode = list(theta = c(0,0), fixed = TRUE)
    )
    
    print(
        all.equal(Sparse(Rst),
                  Sparse(prec(fit3)))
    )

}
