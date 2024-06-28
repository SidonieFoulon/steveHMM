require(steveHMM)
source("../umbrella_example.r")

N <- 500
a <- seq(0,1,length = N)
b <- seq(0,1,length = N)
L <-matrix(0, N, N)
for(i in 1:N) for(j in 1:N) L[i,j] <- neg_log_likelihood( c(a[i],b[j]), X.parapluie, modele.parapluie)

N <- 100
X <- matrix(0, nrow = N, ncol = 2)   
lik <- numeric(N)
for(i in 1:N) {
  set.seed(i)
  theta <- runif(2)
  qnem <- QNEM(theta, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1))
  X[i,] <- as.vector(qnem$theta)
  lik[i] <- qnem$neg.ll
}
# table( round(lik, 2) )

pdf("umbrella_kilehood.pdf")
gaston::lik.contour(a, b, L, xlab = "a", ylab = "b", levels = c(33.33, 33.7, 34.5, 36.6, 38.7, 38.82, 40, 50, 70))
points(X, pch = 4, cex = 2, col = "blue")
dev.off()

