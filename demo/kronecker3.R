library(INLAtools)
R1 <- Matrix(crossprod(diff(diag(4))))
R1
m1 <- cgeneric("generic0", R = R1, param = c(1, NA),
  scale = FALSE, useINLAprecomp = FALSE)
R2 <- Matrix(crossprod(diff(diag(3))))
R2
m2 <- cgeneric("generic0", R = R2, param = c(1, NA),
  scale = FALSE, useINLAprecomp = FALSE)
m3 <- cgeneric("iid", n = 2, param = c(1, 0.5),
  useINLAprecomp = FALSE)
prec(m3, theta = 0.0)
multi123 <- multi_generic_model(
  list(m1 = m1, m2 = m2, m3 = m3),
  useINLAprecomp = FALSE
)
R321 <- Sparse(kronecker(kronecker(Diagonal(2),R2),R1))
R321
all.equal(R321, Sparse(prec(multi123, theta = 0.0)))
if(!is.na(packageCheck("inlabru", "2.13.0.9005"))) {
  print(multi123$mapper)
}
