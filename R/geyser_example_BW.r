## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

p.stationnaire <- function(a,b) c( (1-b), (1-a)*(1-b), a ) / (2 - 2*b + a*b)

modele.geyser <- function(theta, obs) {
  # parametre de la matrice de transition
  a <- theta[1]
  b <- theta[2]
  trans <- matrix(c(0,1-a,a,1,0,0,1-b,0,b), nrow = 3, byrow = TRUE)
  #colnames(trans) <- rownames(trans) <- name.S

  # parametre de l'émission
  c <- theta[3]
  p.emiss <- rbind(ifelse(obs == 1, 1, 0), ifelse(obs == 1, 0, c), ifelse(obs == 1, 0, 1-c))
  #rownames(p.emiss) <- name.S

  # etat stationnaire, ne dépend pas de a
  pi <- p.stationnaire(a,b)

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.geyser <- function(obs, backward) {
  l <- ncol(backward$phi)

  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  a <- (sum(backward$delta[1,3,])) / (l-1)

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  b <- (sum(backward$delta[3,3,])) / (l-1)

  #c = proba de durée >= 3min chez les "Longs" (->emiss)
  c <- sum(backward$phi[2,which(obs==2)]) / l

  c(a,b,c)
}

if(TRUE) {
# Exemple d'utilisation :
#source("forward.r")
#source("backward.r")
#source("EM.r")

# nos observations :
library(MASS)
X <- ifelse(geyser$duration < 3, 1,2)

# nos paramètres a et b d'initialisation :
par <- c(a = 0.3, b = 0.15)

mod <- modele.parapluie(theta = par, obs = X)
fo <- forward(mod)
ba <- backward(mod, fo)

em <- EM(par, X, modele.parapluie, M.step.parapluie, 20)

}

