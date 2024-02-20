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

test.init <- function(obs, it){
  qN <- vector(length = it)
  BW <- vector(length = it)
  SQ <- vector(length = it)

  for(i in 1:it){
    a <- runif(1)
    b <- runif(1)
    par <- c(a = a, b = b)

    f <-function(theta) neg_log_likelihood(theta, obs, modele.parapluie)
    traceCauchy <- captureThetaCauchy(par, f)
    qN[i] <- ncol(traceCauchy) - 1

    em <- EM(par, obs, modele.parapluie, M.step.parapluie, 40)
    em_theta <- em$Theta[, apply(em$Theta, 2, function(x) !all(is.na(x)))]
    BW[i] <- ncol(em_theta)

    squarem <- SQUAREM(par, obs, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 40)
    sq_theta <- squarem[, apply(squarem, 2, function(x) !all(is.na(x)))]
    SQ[i] <- ncol(sq_theta)
  }
 return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

test.time <- function(theta, obs, it){


  #BW
  debut <- Sys.time()
  for(i in 1:it){
    X <- sample(obs)

    em <- EM(theta, X, modele.parapluie, M.step.parapluie, 40)

  }
  fin <- Sys.time()
  BW <- fin - debut

  #SQUAREM

  debut <- Sys.time()
  for(i in 1:it){
    X <- sample(obs)

    squarem <- SQUAREM(par, obs, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 40)
  }
  fin <- Sys.time()
  SQ <- fin - debut

  #qN
  debut <- Sys.time()
  for(i in 1:it){
    X <- sample(obs)

    f <-function(par) neg_log_likelihood(par, X, modele.parapluie)
    traceCauchy <- captureThetaCauchy(theta, f)


  }
  fin <- Sys.time()
  qN <- fin - debut



  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}



if(FALSE) {
# Exemple d'utilisation :

# nos observations :
X <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1,
    2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2)

# nos paramètres a et b d'initialisation :
par <- c(a = 0.3, b = 0.15)

# maximisation directe avec calcul du forward en probabilites conditionnelles
f <-function(theta) neg_log_likelihood(theta, X, modele.parapluie)
optim( par, f, method = "L-B", lower = c(0,0)+0.01, upper = c(1,1))
trace <- captureTheta(par, f)
traceCauchy <- captureThetaCauchy(par, f)

# maximisation directe avec calcul du forward en probabilites jointes
f2 <-function(theta) neg_log_likelihood_Thompson(theta, X, modele.parapluie)
optim( par, f2, method = "L-B", lower = c(0,0)+0.1, upper = c(1,1)-0.1 )
trace2 <- captureTheta(par, f2)
traceCauchy2 <- captureThetaCauchy(par, f2)

# EM
em <- EM(par, X, modele.parapluie, M.step.parapluie, 40)

# SQUAREM
squarem <- SQUAREM(par, X, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 40)

par(mfrow = c(2,1))
plot( em$Theta[1,1:19] , type = "o", ylim = range(em$Theta[1,1:19], squarem[1,1:23], traceCauchy[1,1:19]), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter a")
lines( squarem[1,1:23], type = "o", col = "#c09140")
lines( traceCauchy[1,1:19], type = "o", col = "#c05140")
legend("topright", cex = 1, legend = c("Baum-Welch", "SQUAREM", "quasi-Newton"), col = c("#4080c0", "#c09140", "#c05140"), pch="o")

plot( em$Theta[2,1:19] , type = "o", ylim = range(em$Theta[2,1:19], squarem[2,1:23], traceCauchy[2,1:19]), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter b")
lines( squarem[2,1:23], type = "o", col = "#c09140")
lines( traceCauchy[2,1:19], type = "o", col = "#c05140")
#legend("topright", legend = c("Baum-Welch", "SQUAREM", "quasi-Newton"), col = c("#4080c0", "#c09140", "#c05140"), pch="o")


#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it <- test.init(X, 5000)

#test du temps
set.seed(28)
temps <- test.time(par, obs = X, it = 5000)


}
