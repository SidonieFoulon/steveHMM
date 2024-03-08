require(steveHMM)
source("examples/umbrella_example.r")
em <- EM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, 200)
QNEM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 200)

for(i in 1:100) { 
  set.seed(i); 
  cat("\n***************", i, "************* \n")
  QNEM(runif(2), X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1))
}

set.seed(7)
theta <- runif(2)
QNEM(theta, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1))

N <- 500
a <- seq(0,1,length = N)
b <- seq(0,1,length = N)
L <-matrix(0, N, N)
for(i in 1:N) for(j in 1:N) L[i,j] <- neg_log_likelihood( c(a[i],b[j]), X.parapluie, modele.parapluie)
gaston::lik.contour(a, b, L)



# ----------------------------------------------
source("examples/geyser_dicho.r")
library(MASS)
X.dich <- ifelse(geyser$duration < 3, 1,2)

# nos paramètres a, b et c d'initialisation :
par.dich <- c(a = 0.31, b = 0.46, c = 0.15, d = 0.9, e = 0.9)
em.dich <- EM(par.dich, X.dich, modele.geyser, M.step.geyser, 200)

QNEM(par.dich, X.dich, modele.geyser, M.step.geyser, lower = rep(0,5), upper = rep(1,5))

for(i in 1:100) {
  set.seed(i);
  cat("\n***************", i, "************* \n")
  QNEM(runif(5), X.dich, modele.geyser, M.step.geyser, lower = rep(0,5), upper = rep(1,5))
}

