library(steveHMM)

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


test.init.gd <- function(obs, it){
  qN <- numeric(length = it)
  BW <- numeric(length = it)
  SQ <- numeric(length = it)

  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    par <- runif(5)

    try(qN[i] <- quasi_newton(par, obs, modele.geyser, lower = rep(0.01,5), upper = rep(0.99,5))$counts[1])

    test.em <- function(par){
      em <- EM(par, obs, modele.geyser, M.step.geyser, max.iter = Inf, trace.theta = FALSE)
      return(em$iter)}
    try(BW[i] <- test.em(par))

    test.sq <- function(par){
      squarem <- SQUAREM(par, obs, modele.geyser, M.step.geyser, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = FALSE)
      return(squarem$iter)}
    try(SQ[i] <- test.sq(par))
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

#observations
library(MASS)
X.dich <- ifelse(faithful$eruptions < 3, 1,2)

#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it_gd <- test.init.gd(X.dich, 5000)
saveRDS(nb_it_gd, "/home/sidonie/Bureau/github/steveHMM/examples/nb_it_geyser_dicho.rds")
