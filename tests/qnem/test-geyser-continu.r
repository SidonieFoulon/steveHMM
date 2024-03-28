require(steveHMM)

source("~/GENOSTATS/steveHMM/examples/geyser_continu.r")

X.geyser <- faithful$eruptions
par.geyser <- c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39)


D <- data.frame(theta1 = numeric(100), theta2 = numeric(100), theta3 = numeric(100), theta4 = numeric(100), 
                theta5 = numeric(100), theta6 = numeric(100), theta7 = numeric(100), theta8 = numeric(100), 
       neg.ll = numeric(100), iter = numeric(100), forwards = numeric(100), backwards = numeric(100))
for(i in 1:100) {
  set.seed(i);
  cat("\n***************", i, "************* \n")
  qn <- QNEM( c(runif(2), runif(3, 0, 6), runif(3, 0, 3)), X.geyser, modele.geyser.continu, M.step.geyser.continu, lower = rep(0,8), upper = c(1,1,rep(Inf,6)) )
  # em <- EM( c(runif(2), runif(3, 0, 6), runif(3, 0, 3)), X.geyser, modele.geyser.continu, M.step.geyser.continu, max.iter = Inf)
  # print(round(em$theta, 2))
  # print(em$iter)
  D[i,1:8] <- qn$theta
  D$neg.ll[i] <- qn$neg.ll
  D$iter[i] <- qn$iter
  D$backwards[i] <- qn$backwards
  D$forwards[i] <- qn$forwards
}

set.seed(1)
qn0 <- quasi_newton( c(runif(2), runif(3, 0, 6), runif(3, 0, 3)), X.geyser, modele.geyser.continu, lower = rep(0,8) + 0.01, upper = c(1,1,rep(Inf,6)) - 0.01, trace = TRUE )

