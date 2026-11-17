
library(INLAtools)

####################################################
### this code contain code to
### 1) check cgeneric_wmodel
### 2) example

### auxiliary functions to map theta -> W and W -> theta 
### theta to W 
th2w <- function(th, i) {
    k <- length(th) + 1L
    w <- double(k)
    w[i] <- 1
    w[-i] <- th
    return(w/sqrt(sum(w^2)))
}

th2w(c(-0.5, 2), 1)
th2w(c(-0.5, 2), 2)
th2w(c(-0.5, 2), 3)

### W to theta
w2th <- function(w,i) {
    w[-i]/w[i]
}

all.equal(c(-0.5, 2),
          w2th(th2w(c(-0.5, 2), 1), 1))

#################################################
### 1) check

## generate a random matrix for model Mc
gupp <- sparseMatrix(
    i = c(1, 2, 3, 4, 4, 5),
    j = c(2, 3, 6, 5, 6, 6),
    x = rnorm(6), dims = c(6,6)
); graph1 <- Sparse(abs(gupp + t(gupp))>0)
graph1
(n1 <- nrow(graph1))

### structure matrix R (Laplacian matrix)
###  R_ii = ni, R_ij = -1 if i~j, 0 otherwise
R <- Sparse(
    Diagonal(n1, rowSums(gupp^2) + 1) + gupp + t(gupp))
R

## cgeneric generic0 model: Q = tau * R
K <- 3
Mc <- cgeneric(
    model = "generic0",
    R = R,
    param = c(1, 0.05), ## prior parameters for tau
    constr = FALSE,
    scale = FALSE
)

graph1
graphMc <- cgeneric_graph(Mc) ## query the graph
graphMc
cgeneric_graph(Mc, optimize = TRUE)

all.equal(cgeneric_Q(Mc, theta = 0), ### evaluate Q using cgeneric code
          R)

cgeneric_initial(Mc)

### use the Mc model to define the W-model
###  Q = [ W^{-1} (o) I_n ] bdiag(Q_1,...,Q_K) [W^{-1} (o) I_n]^T
##     = \sum_k W^{-1}[,k](W^{-1}[,k])'(o) Q_k
Wmodel <- cgeneric_wmodel(
    model = Mc,
    K = K,
    lambda = 1
)

Wmodel

K*(K-1) + K
length(cgeneric_initial(Wmodel))

cgeneric_initial(Wmodel)

cgeneric_graph(Wmodel)

cgeneric_mu(Wmodel)

(th1 <- rep(0, K))
(thM <- rnorm(K*(K-1)))

W <- t(sapply(1:K, function(k)
    th2w(thM[(K-1)*(k-1)+1:(K-1)],k)))
W

iW <- solve(W)
iW

## Q from cgeneric
Qth <- cgeneric_Q(Wmodel, theta = c(thM, th1))

## A R funciton to build Q for checking
Qfn <- function(M, lQ) {
    n <- ncol(lQ[[1]])
    MI <- kronecker(M, diag(n))
    MI %*% bdiag(lQ) %*% t(MI)
}

### ordering theta for Q_k
othfn <- function(x) {
    k <- length(x)
    y <- x
    s <- exp(x[k])
    for (i in (k-1):1) {
        s <- s + exp(x[i])
        y[i] <- log(s)
    }
    return(y)
}

or <- othfn(th1)
cbind(th1, or)

## evaluate each Q_k
lQ <- lapply(or, function(x) cgeneric_Q(Mc, theta = x))
## build Q
Qor <- Qfn(iW, lQ)

## compare
all.equal(Sparse(Qth), Sparse(Qor))

##################################################################
### 2) example
rm(list = setdiff(ls(all = TRUE), ### this is against R club rules! but...
                  c("th2w", "w2th", "othfn"))) 
ls()

library(fmesher)
library(INLAtools)
library(INLAspacetime)
library(INLA)

mesh <- fm_mesh_1d(0:299)
cg_s1 <- cgeneric_sspde(
    mesh = mesh, alpha = 2,
    control.priors = list(
        prange = c(3, 0.1), ## P(range < 3) = 0.1
        psigma = c(1, NA) ## fix sigma! 
    )
)

cgeneric_initial(cg_s1)

image(cgeneric_graph(cg_s1))
image(cgeneric_Q(cg_s1, theta = 1:0))

K <- 3
cg_wm <- cgeneric_wmodel(
    model = cg_s1,
    K = K,
    lambda = 2 ## dnorm(x = abs(w_ij), mu = 1/sqrt(K), sigma = lambda)
)

th2w(cgeneric_initial(cg_wm)[1:(K-1)], 1) == (1/sqrt(K))

(th1 <- K:1+1)
lranges <- othfn(th1)
exp(lranges) ## range_k
thM <- rnorm(K * (K-1))

W <- t(sapply(1:K, function(k)
    th2w(thM[(K-1)*(k-1)+1:(K-1)],k)))
W

length(c(thM, th1))

Qth <- cgeneric_Q(cg_wm, theta = c(thM, th1))

image(Qth)

##xx <- inla.qsample(n = 1, Q = Qth)[, 1]
zz <- sapply(lranges, function(lr)
    scale(inla.qsample(n = 1, Q = cgeneric_Q(cg_s1, theta = lr))[,1]))

cov(zz)

xx <- zz %*% t(W)

plot.ts(cbind(zz, xx))

cov(xx)
cor(xx)
tcrossprod(W)

summary(xx)

dat <- data.frame(
    i = 1:cg_wm$f$n,
    y = as.vector(xx)
)

fit <- inla(
    formula = y ~ 0 + f(i, model = cg_wm),
    data = dat,
    control.family = list(
        hyper = list(prec = list(initial = 10, fixed = TRUE))
    )
)

cbind(true = c(thM, th1),
      fit$mode$theta)

Wfit <- t(sapply(1:K, function(k)
    th2w(fit$mode$theta[(K-1)*(k-1)+1:(K-1)],k)))

tcrossprod(W)
tcrossprod(Wfit)
cor(matrix(xx, cg_s1$f$n)) ## data cross-correlation 
