require(steveHMM)
source("~/GENOSTATS/steveHMM/examples/umbrella_example.r")
em <- EM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, 200)
QNEM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 200)

D <- data.frame(theta1 = numeric(100), theta2 = numeric(100), 
     neg.ll = numeric(100), iter = numeric(100), forwards = numeric(100), backwards = numeric(100))
for(i in 1:100) { 
  set.seed(i); 
  cat("\n***************", i, "************* \n")
  qn <- QNEM(runif(2), X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1))
  D[i,1:2] <- qn$theta
  D$neg.ll[i] <- qn$neg.ll
  D$iter[i] <- qn$iter
  D$backwards[i] <- qn$backwards
  D$forwards[i] <- qn$forwards
}

D0 <- readRDS("test-umbrella.rds")
cat("max diff theta1 =", max( abs(D$theta1 - D0$theta1) ), "\n")
cat("max diff theta2 =", max( abs(D$theta2 - D0$theta2) ), "\n")
cat("forwards :", sum(D$forwards), "vs", sum(D0$forwards), "\n")
cat("backwards :", sum(D$backwards), "vs", sum(D0$backwards), "\n")


if(FALSE) {

N <- 500
a <- seq(0,1,length = N)
b <- seq(0,1,length = N)
L <-matrix(0, N, N)
for(i in 1:N) for(j in 1:N) L[i,j] <- neg_log_likelihood( c(a[i],b[j]), X.parapluie, modele.parapluie)
gaston::lik.contour(a, b, L)

}

