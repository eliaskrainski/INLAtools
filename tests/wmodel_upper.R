## Preserving the cross-correlation 
## Z: nxm matrix, Z_{ij} ~ N(0,1),
## W: mxm matrix, cov(ZW') = W'W 
## V(L_k Z.k) = L_k L_k' = V_k
##  X = [L_1 Z.1 : . . . : L_k Z.k ] W'
##  x = vec(X), stacked columns of X
## V(x) = [ W (o) I_n ] bdiag(V_1, ... , V_k) [ W' (o) I_n ]
##      = \sum_k (W.k W.k') (o) V_k
## Q(x) = [ M (o) I_n ] bdiag(Q_1, ... , Q_k) [ M' (o) I_n ]
##      = \sum_k (M.k M.k') (o) Q_k
## with Q_k = V_k^{-1}, M = W^{-1}

library(Matrix)

### 1) check the computation of Q
## play with
iW <- matrix(c(2,-1,-3,4),2)
iW

## A spd matrix
q1 <- matrix(c(4,-3,2, -3,3,-1, 2,-1,2), 3)
q1
chol(q1)

(n1 <- ncol(q1))

## as upper
q1u <- q1*(upper.tri(q1, diag = TRUE))
q1u

lq <- list(q1, q1*3)
lqu <- list(q1u, q1u*3)
qk <- bdiag(lq)
lqk <- bdiag(lqu)
lqk

I1 <- Diagonal(n = n1)
iWI1 <- kronecker(iW, I1)
iWI1

a1 <- iWI1 %*% qk
a1

t(iWI1)
a <- a1 %*% t(iWI1)
a

myQfn <- function(M,lQu) {
### see wmodel.R for checks
### M   : square matrix
### lQu : list of upper matrices 
    K <- ncol(M)
    n <- ncol(lQu[[1]])
    qcompletefn <- function(u) {
        u + t(u) - diag(diag(u))
    }        
    iupp <- upper.tri(diag(n), diag = TRUE)
    out <- matrix(0, n*K, n*K)
    for(i in 1:K) {
        ii <- (i-1)*n + 1:n
        for(j in i:K) {
            jj <- (j-1)*n + 1:n
            if(i==j) {
                for(l in 1:K)
                    out[ii,jj][iupp] <- out[ii,jj][iupp] +
                        M[i,l]*lQu[[l]][iupp]*M[j,l]
            } else {
                for(l in 1:K)
                    out[ii,jj] <- out[ii,jj] + M[i,l]*qcompletefn(lQu[[l]])*M[j,l]
            }
        }
    }
    return(out)
}

au <- as.matrix(a);
au[lower.tri(a)] <- 0
au

all.equal(au, myQfn(iW, lqu))

### 2) the actual W matrix
th2w <- function(th, i) {
    k <- length(th) + 1L
    w <- double(k)
    w[i] <- 1
    w[-i] <- th
    return(w/sqrt(sum(w^2)))
}


### see wmodel.R for (more) checks
rwfn <- function(K)
    t(sapply(1:K, function(i) th2w(rnorm(K-1),i)))
rq2fn <- function(n)
    crossprod(matrix(rnorm(n*n),n))

kk = 5; nn = 10
summary(replicate(20, {
    iW <- solve(rwfn(kk))
    iWI <- kronecker(iW, diag(nn));
    lq <- lapply(1:kk, function(j) rq2fn(nn))
        lqu <- lapply(lq, function(x) x*upper.tri(x,diag=TRUE))
    a <- iWI %*% as.matrix(bdiag(lq)) %*% t(iWI)
    au <- a; au[lower.tri(a)] <- 0
    b <- myQfn(iW, lqu)
    mean((au-b)^2)
}))


### see wmodel.R for more
