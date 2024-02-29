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
    em <- EM(theta, X, modele.parapluie, M.step.parapluie, 1000)
  }
  fin <- Sys.time()
  BW <- fin - debut

  #SQUAREM

  debut <- Sys.time()
  for(i in 1:it){
    X <- sample(obs)
    squarem <- SQUAREM(theta, X, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 1000)
  }
  fin <- Sys.time()
  SQ <- fin - debut

  #qN
  debut <- Sys.time()
  for(i in 1:it){
    X <- sample(obs)
    qN <- capture_quasi_newton(theta, X, modele.parapluie, lower = c(0,0)+0.01, upper = c(1,1)-0.01, trace = TRUE)


  }
  fin <- Sys.time()
  qN <- fin - debut



  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}


# nos observations :
X.parapluie <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1,
    2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2)


if(FALSE) {
# Exemple d'utilisation :

# nos paramètres a et b d'initialisation :
par.parapluie <- c(a = 0.3, b = 0.15)

# maximisation directe avec calcul du forward en probabilites conditionnelles
f <-function(theta) neg_log_likelihood(theta, X.parapluie, modele.parapluie)
optim( par.parapluie, f, method = "L-B", lower = c(0,0)+0.01, upper = c(1,1))
trace <- captureTheta(par.parapluie, f)
traceCauchy <- captureThetaCauchy(par.parapluie, f)
fulltrace <- capture_optim( par.parapluie, f, method = "L-B", lower = c(0,0)+0.01, upper = c(1,1))

# EM
em <- EM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, 200)

# SQUAREM
squarem <- SQUAREM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 200)

#quasi Newton
quasi_newton(par.parapluie, X.parapluie, modele.parapluie, lower = c(0,0)+0.01, upper = c(1,1), trace = TRUE)
qN <- capture_quasi_newton(par.parapluie, X.parapluie, modele.parapluie, lower = c(0,0)+0.01, upper = c(1,1), trace = TRUE)
for(i in 1:ncol(qN)) {
  if( sqrt(sum((qN[,i] - qN[,i-1])^2)) < 1e-5) print(i)
}
qN <- qN[, 1:22]

#plots
par(mfrow = c(2,1))
plot( em$Theta[1,] , type = "o", xlim = c(0, max(ncol(em$Theta), ncol(squarem), ncol(qN))), ylim = range(em$Theta[1,], squarem[1,], qN[1,]), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter a")
abline(v = ncol(em$Theta), col = "#4080c0")
lines( squarem[1,], type = "o", col = "#c09140")
abline(v = ncol(squarem), col = "#c09140")
lines( qN[1,], type = "o", col = "#c05140")
abline(v = ncol(qN), col = "#c05140")
legend("topright", cex = 1, legend = c("Baum-Welch", "SQUAREM", "quasi-Newton"), col = c("#4080c0", "#c09140", "#c05140"), pch="o")

plot( em$Theta[2,] , type = "o", xlim = c(0, max(ncol(em$Theta), ncol(squarem), ncol(qN))), ylim = range(em$Theta[2,], squarem[2,], qN[2,]), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter a")
abline(v = ncol(em$Theta), col = "#4080c0")
lines( squarem[2,], type = "o", col = "#c09140")
abline(v = ncol(squarem), col = "#c09140")
lines( qN[2,], type = "o", col = "#c05140")
abline(v = ncol(qN), col = "#c05140")


#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it <- test.init(X, 5000)

#test du temps
set.seed(28)
temps <- test.time(par.parapluie, obs = X.parapluie, it = 5000)

}
