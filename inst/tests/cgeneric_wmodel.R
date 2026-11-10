
library(INLAtools)

### theta to W map
th2w <- function(th, i) {
    k <- length(th) + 1L
    w <- double(k)
    w[i] <- 1
    w[-i] <- th
    return(w/sqrt(sum(w^2)))
}

### W to theta map
w2th <- function(w,i) {
    w[-i]/w[i]
}

## generate a random matrix for model Mc
gupp <- sparseMatrix(
    i = c(1, 2, 3, 4, 4, 5),
    j = c(2, 3, 6, 5, 6, 6),
    x = rnorm(6), dims = c(6,6)
); graph1 <- Sparse(abs(gupp + t(gupp))>0)
graph1
(n1 <- nrow(graph1))

Q0 <- Sparse(
    Diagonal(n1, rowSums(gupp^2) + 1) + gupp + t(gupp))
Q0

## cgeneric generic0 model
K <- 3
Mc <- cgeneric(
    model = "generic0",
    R = Q0,
    param = c(1, 0.05),
    constr = FALSE,
    scale = FALSE
)

graph1
graphMc <- cgeneric_graph(Mc)
graphMc
cgeneric_graph(Mc, optimize = TRUE)

Q0
cgeneric_Q(Mc, theta = 0)

## use the Mc model to define the W-model
Wmodel <- cgeneric_wmodel(
    model = Mc,
    K = K,
    lambda = 1
)

Wmodel

cgeneric_initial(Wmodel)

cgeneric_graph(Wmodel)

cgeneric_mu(Wmodel)

(th <- rnorm(K*K))

W <- t(sapply(1:K, function(k)
    th2w(th[(K-1)*(k-1)+1:(K-1)],k)))
W
iW <- solve(W)

## Q from cgeneric
Qth <- cgeneric_Q(Wmodel, theta=th)

## build Q function
Qfn <- function(M, lQ) {
    n <- ncol(lQ[[1]])
    MI <- kronecker(M, diag(n))
    MI %*% bdiag(lQ) %*% t(MI)
}

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

or <- othfn(th[K*(K-1)+1:K])
cbind(th[K*(K-1)+1:K], or)

## evuate each Q 
lQ <- lapply(or, function(x) cgeneric_Q(Mc, theta = x))
Qor <- Qfn(iW, lQ)

## compare
all.equal(Sparse(Qth), Sparse(Qor))
