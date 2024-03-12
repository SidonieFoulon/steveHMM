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

test.time.gd <- function(theta, obs, it){
  #qN

  cpt <- 0
  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    X <- sample(obs)
    debut <- Sys.time()
    try(quasi_newton(theta, X, modele.geyser, trace = FALSE, lower = rep(0.01,5), upper = rep(0.99,5)))
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  qN <- cpt

  #BW
  cpt <- 0
  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    X <- sample(obs)
    debut <- Sys.time()
    try(EM(theta, X, modele.geyser, M.step.geyser, max.iter = Inf, trace.theta = FALSE))
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  BW <- cpt

  #SQUAREM

  cpt <- 0
  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    X <- sample(obs)
    debut <- Sys.time()
    try(SQUAREM(theta, X, modele.geyser, M.step.geyser, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = FALSE))
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  SQ <- cpt



  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

#observations
library(MASS)
X.dich <- ifelse(faithful$eruptions < 3, 1,2)

# nos paramètres a, b et c d'initialisation :
par.dich <- c(a = 0.31, b = 0.46, c = 0.15, d = 0.9, e = 0.9)

#test du temps
set.seed(28)
tps_gd <- test.time.gd(par.dich, obs = X.dich, it = 1000)
tps_gd
saveRDS(tps_gd, "/home/sidonie/Bureau/github/steveHMM/examples/tps_geyser_dicho.rds")
