## Preserving the cross-correlation 
## Z: nxm matrix, Z_{ij} ~ N(0,1),
## W: mxm matrix, cov(ZW') = W'W 
## V(L_k Z.k) = L_k L_k' = V_k
##  X = [L_1 Z.1 : . . . : L_k Z.k ] W'
##  x = vec(X), stacked columns of X
## V(x) = [ W (o) I_n ] bdiag(V_1, ... , V_k) [ W' (o) I_n ]
##      = \sum_k (W.k W.k') (o) V_k
## Q(x) = [ M (o) I_n ] bdiag(Q_1, ... , Q_k) [ M' (o) I_n ]
##      = \sum_k ( M.k M.k') (o) Q_k   { ALTERNATIVE }
## with Q_k = V_k^{-1}, M = W^{-1}
## ALTERNATIVE way to compute Q as
## a sum of kroneker products
##    Q(x) = \sum_k ( M[,k] (M[,k])' ) (o) Q_k

library(Matrix)

### 1) check the computation of Q
## play with W^{-1} (some visual easier looking) 
iW <- matrix(c(2,-3,-1,4), 2)
iW

solve(iW) 
rowSums(solve(iW)^2) ## not quite, but just to play with

## (visual easier looking) spd matrix
r1 <- matrix(c(4,-3,2, -3,3,-1, 2,-1,2), 3)
r1
chol(r1)
n1 <- ncol(r1)

## bdiag of r1 and r2*2
lQ <- list(r1, r1*3)
qk <- bdiag(lQ)
qk

I1 <- Diagonal(n = n1)
iWI1 <- kronecker(iW, I1)

iW
iWI1
qk
(a1 <- iWI1 %*% qk)

a1
t(iWI1)
a1 %*% t(iWI1)

Qfn1 <- function(M,lQ) {
    K <- ncol(M)
    n <- ncol(lQ[[1]])
    out <- matrix(0, K*n, K*n)
    for(k in 1:K) {
        wwk <- tcrossprod(M[,k], M[,k])
        out <- out + kronecker(wwk, lQ[[k]])
    }
    return(out)
}

Qfn2 <- function(M,lQ,only1=FALSE) {
    K <- ncol(M)
    n <- ncol(lQ[[1]])
    out <- aux <- matrix(0, n*K, n*K)
### 1st) multiply each Q_k with each element of W
###      |  M[1,1]Q_1  M[1,2]Q_2  ...  M[1,K]Q_K  |
### Q* = |  M[2,1]Q_1  M[2,2]Q_2  ...  M[2,K]Q_K  |
###      |     ...           ...          ...     |
###      |  M[k,1]Q_1  M[k,2]Q_2  ...  M[K,K]Q_K  |
    for(i in 1:K) {
        ii <- (i-1)*n + 1:n
        for(j in 1:K) {
            jj <- (j-1)*n + 1:n
            aux[ii,jj] <- M[i,j] * lQ[[j]]
        }
    }
    if(only1) return(aux)
### 2nd) Q : Q[ii,jj] = sum_i sum_j Q*[ii,jj]M[j,i]
    ## combined 1st and 2nd steps loops
    for(i in 1:K) {
        ii <- (i-1)*n + 1:n
        for(j in 1:K) {
            jj <- (j-1)*n + 1:n
            for(l in 1:K)
                out[ii,jj] <- out[ii,jj] + M[i,l]*lQ[[l]]*M[j,l]
        }
    }
    return(out)
}

all.equal(as.matrix(a1),Qfn2(iW,lQ,TRUE))

b <- as.matrix(a1 %*% t(iWI1))
stopifnot(all.equal(b, Qfn1(iW,lQ)))
stopifnot(all.equal(b, Qfn2(iW,lQ)))

### 2) the actual W matrix
th2w <- function(th, i) {
    k <- length(th) + 1L
    w <- double(k)
    w[i] <- 1
    w[-i] <- th
    return(w/sqrt(sum(w^2)))
}


if(FALSE) {
    
    rW <- replicate(1000, th2w(rnorm(3), 2))
    summary(colSums(rW^2))
    
    par(mfrow = c(2,2))
    for(i in 1:4)
        hist(rW[i,])
    
    rm(rW)

}

if(FALSE) {

    rwfn <- function(K)
        t(sapply(1:K, function(i) th2w(rnorm(K-1),i)))
    summary(replicate(1000, sum(rowSums(rwfn(3)^2))))
    summary(replicate(1000, sum(rowSums(rwfn(4)^2))))
    summary(replicate(1000, sum(rowSums(rwfn(5)^2))))
    summary(replicate(1000, sum(rowSums(rwfn(15)^2))))
    
    r2fn <- function(n)
        crossprod(matrix(rnorm(n*n),n))
    summary(replicate(1000, sum(diag(chol(r2fn(10))))))
    summary(replicate(1000, sum(diag(chol(r2fn(20))))))
    summary(replicate(1000, sum(diag(chol(r2fn(30))))))

    kk = 5; nn = 10
    summary(replicate(20, {
        iW <- solve(rwfn(kk))
        iWI <- kronecker(iW, diag(nn));
        lq <- lapply(1:kk, function(j) r2fn(nn))
        a <- iWI %*% as.matrix(bdiag(lq)) %*% t(iWI)
        b <- Qfn1(iW, lq)
        mean((a-b)^2)
    }))
    
    summary(replicate(20, {
        iW <- solve(rwfn(kk))
        iWI <- kronecker(iW, diag(nn));
        lq <- lapply(1:kk, function(j) r2fn(nn))
        a <- iWI %*% as.matrix(bdiag(lq)) %*% t(iWI)
        b <- Qfn2(iW, lq)
        mean((a-b)^2)
    }))
    
}

###############################################################
###  HAVE A LOOK AT THE MARGINAL COVARIANCE of samples
###############################################################

K <- 3
W <- t(sapply(1:K, function(i) th2w(c(i/K,-i/(i+2)),i)))
W

tcrossprod(W)
iW <- solve(W)
iW

## Model 1 graph
library(INLAtools)
graph1 <- sparseMatrix(
    i = c(1, 2, 3, 4, 4, 5),
    j = c(2, 3, 6, 5, 6, 6), dims = c(6,6)
); graph1 <- Sparse(graph1 + t(graph1))
graph1
(n1 <- nrow(graph1))

Q0 <- Sparse(
    Diagonal(n1, rowSums(graph1) + .1) - graph1)
Q0

cov2cor(chol2inv(chol(as.matrix(Q0))))

### sampling check with W
Z <- matrix(rnorm(K*1e4), 1e4)
tcrossprod(W)
cov(Z %*% t(W))
cov(t(W %*% t(Z)))
cov(t(solve(iW,t(Z))))

### with Q
zz <- matrix(rnorm(n1*1e4), 1e4)
lQ0 <- as.matrix(chol(Q0))
V0 <- chol2inv(lQ0)
lV0 <- chol(V0)

V0
cov(zz %*% lV0)
cov(t(t(lV0) %*% t(zz)))
cov(t(backsolve(lQ0, t(zz))))

### The multivariate
I1 <- Diagonal(n = n1)
iWI1 <- kronecker(iW, I1)
round(iWI1,2)
round(t(iWI1), 2)

## test having each Q_k = Q0
## so the marginal (between rows) are the same
Qlsame <- lapply(1:K, function(k) Q0) 

### cholesky of Q_k (not of V_k)
LQlsame <- lapply(Qlsame, chol)
lQlsame <- lapply(LQlsame, as.matrix)

### drawn big number of samples to approximate the theoretical
nsim <- 10000
xx <- lapply(1:nsim, function(i) {
    x <- matrix(rnorm(K*n1),n1)
    for(k in 1:K)
        x[,k] <- backsolve(lQlsame[[k]], x[,k])
    x %*% W
})

tcrossprod(W)
cor(do.call('rbind', xx))

chol2inv(chol(as.matrix(Q0)))
cov(t(do.call('cbind', xx)))

## different Q: now the between-row marginals are a mixture
## (preserve the "pattern" as only the variance of Q_k differ)
Qlk <- lapply(1:K, function(k) Q0 * k) 
Qk <- Sparse(bdiag(Qlk))

round(iWI1, 2)
Qk
round(iWI1 %*% Qk, 2)

Q <- Sparse(iWI1 %*% Qk %*% t(iWI1))
Q
LQ <- chol(Q)
lQ <- as.matrix(LQ)

xxQ <- t(backsolve(lQ, matrix(rnorm(n1*K*nsim),n1*K)))

library(fields)

fcol <- tim.colors(40)
cbk <- seq(-1,1,length=41)

par(mfrow = c(2,2), mar = c(.5,.5,.5,5), mgp = c(1.5,0.5,0))
image.plot(chol2inv(lQ), axes = FALSE)
image.plot(cov(xxQ), axes = FALSE)
image.plot(cov2cor(chol2inv(lQ)), axes = FALSE, breaks=cbk, col=fcol)
image.plot(cor(xxQ), axes = FALSE, breaks=cbk, col=fcol)

tcrossprod(W)
chol2inv(chol(as.matrix(Q0)))
cov2cor(chol2inv(chol(as.matrix(Q0))))
