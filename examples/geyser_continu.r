## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

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

  # etat stationnaire
  pi <- c( c = (1-b), l = (1-a)*(1-b), ls = a ) / (2 - 2*b + a*b)

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.geyser.continu <- function(obs, backward) {
  l <- ncol(backward$phi)

  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  a <- sum(backward$delta[1,3,-1]) / sum(backward$phi[1,-l])

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  b <- sum(backward$delta[3,3,-1]) / sum(backward$phi[3,-l])

  #esperance de la loi normale pour les etats "court" et "Long", la proba pour l'état "Long stable" en découlera (->emiss)
  muc <- sum(backward$phi[1,] * obs) / sum(backward$phi[1,])

  mul <- sum(backward$phi[2,] * obs) / sum(backward$phi[2,])

  muls <- sum(backward$phi[3,] * obs) / sum(backward$phi[3,])

  #ecart-type de la loi normale pour les etats "court" et "Long", la proba pour l'état "Long stable" en découlera (->emiss)
  varc <- (sum(backward$phi[1,] * (obs**2)) / sum(backward$phi[1,])) - (muc**2)
  sdc <- sqrt(ifelse(varc > 0, varc, 0))

  varl <- (sum(backward$phi[2,] * (obs**2)) / sum(backward$phi[2,])) - (mul**2)
  sdl <- sqrt(ifelse(varl > 0, varl, 0))

  varls <- (sum(backward$phi[3,] * (obs**2)) / sum(backward$phi[3,])) - (muls**2)
  sdls <- sqrt(ifelse(varls > 0, varls, 0))

  r <- c(a = a, b = b, muc = muc, mul = mul, muls = muls, sdc = sdc, sdl = sdl, sdls = sdls)
  r[is.na(r)] <- 0
  r
}

if(FALSE) {
# nos observations continues :
library(MASS)
X <- geyser$duration

# nos paramètres a et b et paramètres des lois normales d'initialisation :
par <- c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39)


# maximisation directe avec calcul du forward en probabilites conditionnelles
f <-function(theta) neg_log_likelihood(theta, X, modele.geyser.continu)
optim( par, f, method = "L-B", lower = c(0,0), upper = c(1,1) )
traceCauchy <- captureThetaCauchy(par, f, low = c(0,0), up = c(1,1))

# EM
em.dich <- EM(par.dich, X.dich, modele.geyser, M.step.geyser, 20)

# SQUAREM
squarem.dich <- SQUAREM(par.dich, X.dich, modele.geyser, M.step.geyser, lower = c(0,0), upper = c(1,1), 20)

# quasi newton avec dérivées
quasi_newton(par, X, modele.geyser.continu, lower = rep(0,5)+1e-5, upper = c(1,1, rep(Inf,3)), trace = TRUE)


}
