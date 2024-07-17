## on note : S in {Dry : 1 ; Rainy : 2}
##           X in {No umbrella : 1 ; Umbrella : 2}

modele.parapluie <- function(theta, obs, name.S = c("Dry", "Rainy")) {
  # parametre de la matrice de transition
  a <- theta[1]
  trans <- matrix( c(1-a, a, a, 1-a), byrow = TRUE, nrow = 2, ncol = 2)
  colnames(trans) <- rownames(trans) <- name.S
  # etat stationnaire, ici ne dépend pas de a ou b
  stat <- c(0.5, 0.5)
  # parametre de l'émission
  b <- theta[2]
  p.emiss <- rbind(ifelse(obs == 1, 1-b, b), ifelse(obs == 1, b, 1-b))
  rownames(p.emiss) <- name.S
  list(trans = trans, pi = stat, p.emiss = p.emiss)
}

M.step.parapluie <- function(obs, backward) {
  l <- ncol(backward$phi)
  # a = proba de changement d'état caché (-> trans)
  a <- (sum(backward$delta[1,2,] + backward$delta[2,1,]))/(l-1)
  # b = proba d'erreur dans le choix No Umbrella/Umbrella (-> emiss)
  b <- (sum(backward$phi[2,which(obs == 1)]) + sum(backward$phi[1,which(obs == 2)]))/l
  c(a,b)
}


# nos observations :
X.parapluie <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1,
    2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2)



# Exemple d'utilisation :

# nos paramètres a et b d'initialisation :
par.parapluie <- c(a = 0.3, b = 0.15)

# EM
em <- EM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, max.iter = 200, trace.theta = TRUE)

# SQUAREM
squarem <- SQUAREM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1),max.iter =  200, trace.theta = TRUE)

#quasi Newton
quasi_newton(par.parapluie, X.parapluie, modele.parapluie, lower = c(0,0)+0.01, upper = c(1,1), trace = TRUE)



