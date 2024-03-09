require(steveHMM)

source("~/GENOSTATS/steveHMM/examples/geyser_dicho.r")
library(MASS)
X.dich <- ifelse(geyser$duration < 3, 1,2)

# nos paramètres a, b et c d'initialisation :
par.dich <- c(a = 0.31, b = 0.46, c = 0.15, d = 0.9, e = 0.9)
em.dich <- EM(par.dich, X.dich, modele.geyser, M.step.geyser, 200)

QNEM(par.dich, X.dich, modele.geyser, M.step.geyser, lower = rep(0,5), upper = rep(1,5))

D <- data.frame(theta1 = numeric(100), theta2 = numeric(100), theta3 = numeric(100), theta4 = numeric(100), theta5 = numeric(100),
       neg.ll = numeric(100), iter = numeric(100), forwards = numeric(100), backwards = numeric(100))
for(i in 1:100) {
  set.seed(i);
  cat("\n***************", i, "************* \n")
  qn <- QNEM(runif(5), X.dich, modele.geyser, M.step.geyser, lower = rep(0,5), upper = rep(1,5))
  D[i,1:5] <- qn$theta
  D$neg.ll[i] <- qn$neg.ll
  D$iter[i] <- qn$iter
  D$backwards[i] <- qn$backwards
  D$forwards[i] <- qn$forwards
}

D0 <- readRDS("test-geyser.rds")
cat("max diff theta1 =", max( abs(D$theta1 - D0$theta1) ), "\n")
cat("max diff theta2 =", max( abs(D$theta2 - D0$theta2) ), "\n")
cat("max diff theta3 =", max( abs(D$theta3 - D0$theta3) ), "\n")
cat("max diff theta4 =", max( abs(D$theta4 - D0$theta4) ), "\n")
cat("max diff theta5 =", max( abs(D$theta5 - D0$theta5) ), "\n")
cat("forwards :", sum(D$forwards), "vs", sum(D0$forwards), "\n")
cat("backwards :", sum(D$backwards), "vs", sum(D0$backwards), "\n")


