library(steveHMM)

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

test.init <- function(obs, it){
  qN <- vector(length = it)
  BW <- vector(length = it)
  SQ <- vector(length = it)

  for(i in 1:it){
    print(i)
    par <- runif(2)

    f <-function(theta) neg_log_likelihood(theta, obs, modele.parapluie)
    qN[i] <- quasi_newton(par, obs, modele.parapluie, lower = c(0,0)+1e-2, upper = c(1,1)-1e-2)$counts[1]

    em <- EM(par, obs, modele.parapluie, M.step.parapluie, max.iter = Inf, trace.theta = FALSE)
    BW[i] <- em$iter

    squarem <- SQUAREM(par, obs, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), max.iter = Inf, trace.theta = FALSE)
    SQ[i] <- squarem$iter
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

# nos observations :
X.parapluie <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1,
                 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2)

#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it <- test.init(X.parapluie, 5000)
saveRDS(nb_it, "/home/sidonie/Bureau/github/steveHMM/examples/nb_it_umbrella.rds")
