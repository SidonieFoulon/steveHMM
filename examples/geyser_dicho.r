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

  c(a, b, c, d, e)
}

test.init <- function(obs, it, max.iter = 1500){
  qN <- vector(length = it)
  BW <- vector(length = it)
  SQ <- vector(length = it)

  for(i in 1:it){
    print(i)
    par <- runif(5, 0, 1)
    print(par)

    print("Debut qN")
    f <-function(theta) neg_log_likelihood(theta, obs, modele.geyser)
    traceCauchy <- captureThetaCauchy(par, f, low = c(0,0,0,0,0) + 0.01, up = c(1,1,1,1,1) - 0.01)
    qN[i] <- ncol(traceCauchy)

    print("Debut EM")
    em <- EM(par, obs, modele.geyser, M.step.geyser, max.iter)
    em_theta <- em$Theta[, apply(em$Theta, 2, function(x) !all(is.na(x)))]
    BW[i] <- ncol(em_theta)

    print("Debut SQ")
    squarem <- SQUAREM(par, obs, modele.geyser, M.step.geyser, lower = c(0,0,0,0,0), upper = c(1,1,1,1,1), max.iter)
    sq_theta <- squarem[, apply(squarem, 2, function(x) !all(is.na(x)))]
    SQ[i] <- ncol(sq_theta)
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

if(FALSE) {
# Exemple d'utilisation :
#source("forward.r")
#source("backward.r")
#source("EM.r")

## Version dichotomisée

# nos observations dichotomisées :
library(MASS)
X.dich <- ifelse(geyser$duration < 3, 1,2)

# nos paramètres a, b et c d'initialisation :
par.dich <- c(a = 0.31, b = 0.46, c = 0.15, d = 0.9, e = 0.9)


# maximisation directe avec calcul du forward en probabilites conditionnelles
f <-function(theta) neg_log_likelihood(theta, X.dich, modele.geyser)
optim( par.dich, f, method = "L-B", lower = c(0,0,0,0,0) + 0.01, upper = c(1,1,1,1,1) - 0.01)
traceCauchy.dich <- captureThetaCauchy(par.dich, f, low = c(0,0,0,0,0)+ 0.01, up = c(1,1,1,1,1) - 0.01)

# maximisation directe avec calcul du forward en probabilites jointes
f2 <-function(theta) neg_log_likelihood_Thompson(theta, X.dich, modele.geyser)
optim( par.dich, f2, method = "L-B", lower = c(0,0,0)+0.01, upper = c(1,1,1)-0.01 )
traceCauchy2.dich <- captureThetaCauchy(par.dich, f2, low = c(0,0,0)+ 0.01, up = c(1,1,1) - 0.01)

# EM
em.dich <- EM(par.dich, X.dich, modele.geyser, M.step.geyser, 200)

# SQUAREM
squarem.dich <- SQUAREM(par.dich, X.dich, modele.geyser, M.step.geyser, lower = c(0,0,0,0,0), upper = c(1,1,1,1,1), 200)


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

#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it.g <- test.init(X.dich, 10)


}

