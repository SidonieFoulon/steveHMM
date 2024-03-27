## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

p.stationnaire <- function(a,b) {
  c( c = (1-b), l = (1-a)*(1-b), ls = a ) / (2 - 2*b + a*b)
}

modele.geyser <- function(theta, obs, name.S = c("Court", "Long", "Long stable")) {
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
  pi_c <- theta[6]
  pi_l <- theta[7]
  pi_ls <- theta[8]
  pi <- c(c = pi_c, l = pi_l, ls = pi_ls)

  list(trans = trans, pi = pi, p.emiss = p.emiss)
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

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.geyser <- function(obs, backward) {
  l <- ncol(backward$phi)
  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  a <- sum(backward$delta[1,3,-1]) / sum(backward$phi[1,-l])

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  b <- sum(backward$delta[3,3,-1]) / sum(backward$phi[3,-l])

  #c = proba de durée >= 3min chez les "Courts" (->emiss)
  c <- sum(backward$phi[1,which(obs==2)]) / sum(backward$phi[1,])

  #d = proba de durée >= 3min chez les "Longs" (->emiss)
  d <- sum(backward$phi[2,which(obs==2)]) / sum(backward$phi[2,])

  #e = proba de durée >= 3min chez les "Longs stables" (->emiss)
  e <- sum(backward$phi[3,which(obs==2)]) / sum(backward$phi[3,])

  #etats initiaux
  pi_c <- backward$phi[1,1]
  pi_l <- backward$phi[2,1]
  pi_ls <- backward$phi[3,1]

  c(a, b, c, d, e, pi_c, pi_l, pi_ls)
}

M.step.geyser.stat <- function(obs, backward) {
  l <- ncol(backward$phi)
  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  a <- sum(backward$delta[1,3,-1]) / sum(backward$phi[1,-l])

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  b <- sum(backward$delta[3,3,-1]) / sum(backward$phi[3,-l])

  #c = proba de durée >= 3min chez les "Courts" (->emiss)
  c <- sum(backward$phi[1,which(obs==2)]) / sum(backward$phi[1,])

  #d = proba de durée >= 3min chez les "Longs" (->emiss)
  d <- sum(backward$phi[2,which(obs==2)]) / sum(backward$phi[2,])

  #e = proba de durée >= 3min chez les "Longs stables" (->emiss)
  e <- sum(backward$phi[3,which(obs==2)]) / sum(backward$phi[3,])

  c(a, b, c, d, e)
}


if(FALSE) {
# Exemple d'utilisation :
#source("forward.r")
#source("backward.r")
#source("EM.r")

## Version dichotomisée

# nos observations dichotomisées :
library(MASS)
X.dich <- ifelse(faithful$eruptions < 3, 1,2)

# nos paramètres a, b et c d'initialisation :
pi.dich <- p.stationnaire(0.31,0.46)
par.dich <- c(a = 0.31, b = 0.46, c = 0.15, d = 0.9, e = 0.9, pi_c = pi.dich[1], pi_l = pi.dich[2], pi_ls = pi.dich[3])


# qN
qn.dich <- quasi_newton(par.dich, X.dich, modele.geyser, lower = rep(0,8) + 0.01, upper = rep(1,8) - 0.01)
qn.dich.trace <- capture_quasi_newton(par.dich, X.dich, modele.geyser, lower = rep(0,8) + 0.01, upper = rep(1,8) - 0.01)

# EM
em.dich <- EM(par.dich, X.dich, modele.geyser, M.step.geyser, max.iter = 1000, trace.theta = TRUE)

# SQUAREM
squarem.dich <- SQUAREM(par.dich, X.dich, modele.geyser, M.step.geyser, lower = rep(0,8), upper = rep(1,8), max.iter = 1000, trace.theta = TRUE)

# QNEM
qnem.dich <- QNEM(par.dich, X.dich, modele.geyser, M.step.geyser, max.iter = 1000, lower = rep(0,8), upper = rep(1,8), trace.theta = TRUE, verbose = TRUE)

#plot
par(mfrow = c(5,1))
em.dich$Theta <- em.dich$Theta[, apply(em.dich$Theta, 2, function(x) !all(is.na(x)))]
squarem.dich <- squarem.dich[, apply(squarem.dich, 2, function(x) !all(is.na(x)))]
plot( em.dich$Theta[1,] , type = "o", ylim = c(0,1), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter a")
lines( squarem.dich[1,], type = "o", col = "#c09140")
lines( traceCauchy.dich[1,], type = "o", col = "#c05140")
legend("topright", cex = 1, legend = c("Baum-Welch", "SQUAREM", "quasi-Newton"), col = c("#4080c0", "#c09140", "#c05140"), pch="o")

plot( em.dich$Theta[2,] , type = "o", ylim = c(0,1), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter b")
lines( squarem.dich[2,], type = "o", col = "#c09140")
lines( traceCauchy.dich[2,], type = "o", col = "#c05140")

plot( em.dich$Theta[3,] , type = "o", ylim = c(0,1), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter c")
lines( squarem.dich[3,], type = "o", col = "#c09140")
lines( traceCauchy.dich[3,], type = "o", col = "#c05140")

plot( em.dich$Theta[4,] , type = "o", ylim = c(0,1), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter d")
lines( squarem.dich[4,], type = "o", col = "#c09140")
lines( traceCauchy.dich[4,], type = "o", col = "#c05140")

plot( em.dich$Theta[5,] , type = "o", ylim = c(0,1), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter e")
lines( squarem.dich[5,], type = "o", col = "#c09140")
lines( traceCauchy.dich[5,], type = "o", col = "#c05140")


}

