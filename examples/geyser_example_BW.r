## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

p.stationnaire <- function(a,b) c( c = (1-b), l = (1-a)*(1-b), ls = a ) / (2 - 2*b + a*b)

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

modele.geyser.continu <- function(theta = c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39), obs, name.S = c("Court", "Long", "Long stable")) { #theta contient a, b + mu, sigma
  # parametre de la matrice de transition
  a <- theta[1]
  b <- theta[2]
  trans <- matrix(c(0,1-a,a,1,0,0,1-b,0,b), nrow = 3, byrow = TRUE)
  colnames(trans) <- rownames(trans) <- name.S

  # parametre de l'émission
  muc <- theta[3]
  mul <- theta[4]
  muls <- theta[5]
  sdc <- theta[6]
  sdl <- theta[7]
  sdls <- theta[8]

  #proba de l'observation par la loi normale
  p.emiss.c <- dnorm(obs, mean = muc, sd = sdc)
  p.emiss.l <- dnorm(obs, mean = mul, sd = sdl)
  p.emiss.ls <- dnorm(obs, mean = muls, sd = sdls)

  p.emiss <- rbind(p.emiss.c, p.emiss.l, p.emiss.ls)
  rownames(p.emiss) <- name.S

  # etat stationnaire, ne dépend pas de a
  pi <- p.stationnaire(a,b)

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.geyser <- function(obs, backward) {
  l <- ncol(backward$phi)
  div0 <- function(x,y) ifelse(y == 0, 0, x / y)
  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  #a <- sum( div0(backward$delta[1,3,-1] , backward$phi[1,-l])) / sum( backward$phi[1,-l])
  a <- div0(sum(backward$delta[1,3,-1]) , sum(backward$phi[1,-l]))

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  #b <- sum( div0(backward$delta[3,3,-1], backward$phi[3,-l])) / sum(backward$phi[3,-l])
  b <- div0(sum(backward$delta[3,3,-1]), sum(backward$phi[3,-l]))


  #c = proba de durée >= 3min chez les "Courts" (->emiss)
  c <- sum(backward$phi[1,which(obs==2)]) / sum(backward$phi[1,])

  #d = proba de durée >= 3min chez les "Longs" (->emiss)
  d <- sum(backward$phi[2,which(obs==2)]) / sum(backward$phi[2,])

  #e = proba de durée >= 3min chez les "Longs stables" (->emiss)
  e <- sum(backward$phi[3,which(obs==2)]) / sum(backward$phi[3,])

  c(a, b, c, d, e)
}

M.step.geyser.continu <- function(obs, backward) {
  l <- ncol(backward$phi)

  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  a <- (sum(backward$delta[1,3,])) / (l-1)

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  b <- (sum(backward$delta[3,3,])) / (l-1)

  #esperance de la loi normale pour les etats "court" et "Long", la proba pour l'état "Long stable" en découlera (->emiss)
  muc <- sum(backward$phi[1,] * obs) / sum(backward$phi[1,])

  mul <- sum(backward$phi[2,] * obs) / sum(backward$phi[2,])

  muls <- sum(backward$phi[3,] * obs) / sum(backward$phi[3,])

  #ecart-type de la loi normale pour les etats "court" et "Long", la proba pour l'état "Long stable" en découlera (->emiss)
  varc <- (sum(backward$phi[1,] * (obs**2)) / sum(backward$phi[1,])) - (muc**2)
  sdc <- sqrt(varc)

  varl <- (sum(backward$phi[2,] * (obs**2)) / sum(backward$phi[2,])) - (mul**2)
  sdl <- sqrt(varl)

  varls <- (sum(backward$phi[3,] * (obs**2)) / sum(backward$phi[3,])) - (muls**2)
  sdls <- sqrt(varls)

  c(a = a, b = b, muc = muc, mul = mul, muls = muls, sdc = sdc, sdl = sdl, sdls = sdls)
}


test.init <- function(obs, it){
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
    em <- EM(par, obs, modele.geyser, M.step.geyser, 100)
    em_theta <- em$Theta[, apply(em$Theta, 2, function(x) !all(is.na(x)))]
    BW[i] <- ncol(em_theta)

    print("Debut SQ")
    squarem <- SQUAREM(par, obs, modele.geyser, M.step.geyser, lower = c(0,0,0,0,0), upper = c(1,1,1,1,1), 100)
    sq_theta <- squarem[, apply(squarem, 2, function(x) !all(is.na(x)))]
    SQ[i] <- ncol(sq_theta)
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

if(TRUE) {
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

## Version continue

# nos observations continues :
library(MASS)
X <- geyser$duration

# nos paramètres a et b et paramètres des lois normales d'initialisation :
par <- c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39)


# maximisation directe avec calcul du forward en probabilites conditionnelles
f <-function(theta) neg_log_likelihood(theta, X, modele.geyser.continu)
optim( par, f, method = "L-B", lower = c(0,0), upper = c(1,1) )
traceCauchy <- captureThetaCauchy(par, f, low = c(0,0), up = c(1,1))

# maximisation directe avec calcul du forward en probabilites jointes
f2 <-function(theta) neg_log_likelihood_Thompson(theta, X, modele.geyser.continu)
optim( par, f2, method = "L-B", lower = c(0,0) + 0.1, upper = c(1,1) - 0.1)
traceCauchy2.dich <- captureThetaCauchy(par.dich, f2, low = c(0,0)+ 0.01, up = c(1,1) - 0.01)

# EM
em.dich <- EM(par.dich, X.dich, modele.geyser, M.step.geyser, 20)

# SQUAREM
squarem.dich <- SQUAREM(par.dich, X.dich, modele.geyser, M.step.geyser, lower = c(0,0), upper = c(1,1), 20)



}

