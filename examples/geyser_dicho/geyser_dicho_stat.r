## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

p.stationnaire <- function(a,b) {
  c( c = (1-b), l = (1-a)*(1-b), ls = a ) / (2 - 2*b + a*b)
}

modele.geyser.stat <- function(theta, obs, name.S = c("Court", "Long", "Long stable")) {
  # parametre de la matrice de transition
  a <- theta[1]
  b <- theta[2]
  trans <- matrix(c(0,1-a,a,1,0,0,1-b,0,b), nrow = 3, byrow = TRUE)
  colnames(trans) <- rownames(trans) <- name.S

  # parametre de l'émission
  c <- theta[3]
  d <- theta[4]
  e <- theta[5]
  p.emiss <- rbind(ifelse(obs == 1, 1-c, c), ifelse(obs == 1, 1-d, d), ifelse(obs == 1, 1-e, e))
  rownames(p.emiss) <- name.S

  # etat stationnaire, ne dépend pas de a
  pi <- p.stationnaire(a,b)

  list(trans = trans, pi = pi, p.emiss = p.emiss, theta = theta)
}

M.step.geyser.stat <- function(obs, backward) {
  l <- ncol(backward$phi)
  
  a0 <- backward$theta[1]
  b0 <- backward$theta[2]
  D12 <- sum(backward$delta[1,2,-1]) + backward$phi[2,1]
  D13 <- sum(backward$delta[1,3,-1]) + backward$phi[3,1]
  D31 <- sum(backward$delta[3,1,-1]) + backward$phi[1,1] + backward$phi[2,1]
  D33 <- sum(backward$delta[3,3,-1])

  repeat {
    if(D12 + D13 > 0)
      a <- (D13 - a0*(1-a0)*b0/(2 - 2*b0 + a0*b0)) / (D12 + D13)
    else
      a <- 1 # pas d'état court

    if(D31 + D33 > 0)
      b <- (D33 - b0*(1-b0)*(a0-2)/(2 - 2*b0 + a0*b0)) / (D31 + D33)
    else
      b <- 0 # pas d'état "long stable"

    if(a >= 0 & a <= 1 & b >= 0 & b <= 1) break
    a0 <- a
    b0 <- b
  }

  #c = proba de durée >= 3min chez les "Courts" (->emiss)
  c <- sum(backward$phi[1,which(obs==2)]) / sum(backward$phi[1,])

  #d = proba de durée >= 3min chez les "Longs" (->emiss)
  d <- sum(backward$phi[2,which(obs==2)]) / sum(backward$phi[2,])

  #e = proba de durée >= 3min chez les "Longs stables" (->emiss)
  e <- sum(backward$phi[3,which(obs==2)]) / sum(backward$phi[3,])

  c(a, b, c, d, e)
}

# nos observations dichotomisées :
X.dich <- ifelse(faithful$eruptions < 3, 1,2)

# nos paramètres a, b et c d'initialisation :
par.dich.stat <- c(a = 0.31, b = 0.46, c = 0.15, d = 0.9, e = 0.9)


if(FALSE) {
library(steveHMM)

# qN
qn.dich <- quasi_newton(par.dich.stat, X.dich, modele.geyser.stat, lower = rep(0,7) + 0.01, upper = rep(1,7) - 0.01)
qn.dich.trace <- capture_quasi_newton(par.dich.stat, X.dich, modele.geyser.stat, lower = rep(0,7) + 0.01, upper = rep(1,7) - 0.01)

# EM
em.dich <- EM(par.dich.stat, X.dich, modele.geyser.stat, M.step.geyser.stat, max.iter = 1000, trace.theta = TRUE)

# SQUAREM
squarem.dich <- SQUAREM(par.dich.stat, X.dich, modele.geyser.stat, M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = 1000, trace.theta = TRUE)

# QNEM
qnem.dich <- QNEM(par.dich.stat, X.dich, modele.geyser.stat, M.step.geyser.stat, max.iter = 1000, lower = rep(0,5), upper = rep(1,5), trace.theta = TRUE, verbose = TRUE)
}

