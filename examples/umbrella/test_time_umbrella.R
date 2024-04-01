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

test.time <- function(theta, obs, it){


  #BW
  cpt <- 0
  for(i in 1:it){
    print(i)
    X <- sample(obs)
    debut <- Sys.time()
    EM(theta, X, modele.parapluie, M.step.parapluie, max.iter = 1000, trace.theta = FALSE)
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  BW <- cpt

  #SQUAREM

  cpt <- 0
  for(i in 1:it){
    print(i)
    X <- sample(obs)
    debut <- Sys.time()
    SQUAREM(theta, X, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), max.iter = 1000, trace.theta = FALSE)
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  SQ <- cpt

  #qN
  cpt <- 0
  for(i in 1:it){
    print(i)
    X <- sample(obs)
    debut <- Sys.time()
    quasi_newton(theta, X, modele.parapluie, trace = FALSE, lower = c(0,0)+1e-2, upper = c(1,1)-1e-2)
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  qN <- cpt



  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

# nos observations :
X.parapluie <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1,
                 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2)


# nos paramètres a et b d'initialisation :
par.parapluie <- c(a = 0.3, b = 0.15)

#test du temps
set.seed(28)
temps <- test.time(par.parapluie, obs = X.parapluie, it = 5000)
temps
