
## Define a model with precision is Q = Q1 (x) Q2
## Use two ways to define it and fit using inla
## We set m1, m2 as the 'generic0',
##  with Q1 = \theta * R and Q2  = 1 * R2

library(INLA)
library(INLAtools)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1 graph
graph <- sparseMatrix(
    i = c(2, 3, 1, 4, 1, 4, 2, 3),
    j = c(1, 1, 2, 2, 3, 3, 4, 4)
)
graph 
(n <- nrow(graph))

R1 <- as(inla.as.sparse(
    Diagonal(n, rowSums(graph)) - graph),
    'symmetricMatrix')
R1[1:min(10, n), 1:min(20, n)]

## model 1 definition
m1 <- cgeneric(
    model = 'generic0',
    R = R1,
    scale = FALSE,
    constr = FALSE,
    param = c(1, 0.5),
    debug = FALSE)

## model 2 graph 
n2 <- 4
R2 <- INLA:::inla.rw1(n2)
R2[1:min(5, n2), 1:min(20, n2)]

## model 2 cgeneric definition
m2 <- cgeneric(
    model = 'generic0',
    R = R2,
    scale = FALSE,
    constr = FALSE,
    param = c(1, 0.0),
    debug = FALSE)

Q21 <- as(inla.as.sparse(kronecker(R2, R1)), 
          'symmetricMatrix')
Q21[1:min(5, n*n2), 1:min(10, n*n2)]

(theta.fixed <- c(0))

cmode <- list(
    theta = theta.fixed,
    fixed = TRUE,
    restart = FALSE)
ccomp <- list(config = TRUE)

dataf <- data.frame(
    y = rep(NA, n*n2),
    i = rep(1:n, each = n2),
    j = rep(1:n2, n))

ires1 <- inla(
    y ~ 0 + f(i, model = "generic0", Cmatrix = R1, group = j, 
              control.group = list(model = 'rw1', scale.model = FALSE)),
    family = "poisson",
    data = dataf,
    control.mode = cmode,
    control.compute = ccomp
)

Qinla1 <- ires1$misc$configs$config[[1]]$Q

all.equal(Sparse(upperPadding(Q21)),
          Sparse(Qinla1))

## kronecker model
k21 <- kronecker(m2, m1)

q21 <- cgeneric_Q(k21, theta = c(theta.fixed))

all.equal(q21, Sparse(Q21))

all.equal(Sparse(Qinla1),
          Sparse(upperPadding(q21)))

ires2 <- inla(
    y ~ 0 + f(i, model = k21),
    data = data.frame(y = NA, i = 1:(n * n2)),
    family = "poisson",
    control.compute = list(config = TRUE),
    control.mode =
        list(theta = c(theta.fixed),
             fixed = TRUE)
)

Qinla2 <- ires2$misc$configs$config[[1]]$Q
all.equal(Sparse(upperPadding(q21)),
          Sparse(Qinla2))

detach("package:INLAtools", unload = TRUE)
library(INLAtools)
