## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

modele.geyser.continu.stat <- function(theta = c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39), obs, name.S = c("Court", "Long", "Long stable")) { #theta contient a, b + mu, sigma
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

  list(trans = trans, pi = pi, p.emiss = p.emiss, theta = theta)
}

M.step.geyser.continu.stat <- function(obs, backward) {
  l <- ncol(backward$phi)

  a0 <- backward$theta[1]
  b0 <- backward$theta[2]
  D12 <- sum(backward$delta[1,2,-1]) + backward$phi[2,1]
  D13 <- sum(backward$delta[1,3,-1]) + backward$phi[3,1]
  D31 <- sum(backward$delta[3,1,-1]) + backward$phi[1,1] + backward$phi[2,1]
  D33 <- sum(backward$delta[3,3,-1]) 
 
  if(D12 + D13 > 0) 
    a <- (D13 - a0*(1-a0)/(2 - 2*b0 + a0*b0)) / (D12 + D13)
  else 
    a <- 1 # pas d'état court

  if(D31 + D33 > 0) 
    b <- (D33 - b0*(1-b0)/(2 - 2*b0 + a0*b0)) / (D31 + D33)
  else 
    b <- 0 # pas d'état "long stable"


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
  r
}


# données geyser
X.geyser <- faithful$eruptions

# nos paramètres a et b et paramètres des lois normales d'initialisation :
par.geyser.stat <- c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39)


if(FALSE) {
library(steveHMM)
# qN
qn.cont <- quasi_newton(par.geyser.stat, X.geyser, modele.geyser.continu.stat, lower = c(.01,.01, rep(0, 6), .01, .01), upper = c(.99,.99, rep(Inf, 6), .99,.99))
qn.cont.trace <- capture_quasi_newton(par.geyser.stat, X.geyser, modele.geyser.continu.stat, lower = c(.01,.01, rep(0, 6), .01, .01), upper = c(.99,.99, rep(Inf, 6), .99,.99))

# EM
em.cont <- EM(par.geyser.stat, X.geyser, modele.geyser.continu.stat, M.step.geyser.continu.stat, trace.theta = TRUE)

# SQUAREM
squarem.cont <- SQUAREM(par.geyser.stat, X.geyser, modele.geyser.continu.stat, M.step.geyser.continu.stat, trace.theta = TRUE)

# QNEM
qnem.cont <- QNEM(par.geyser.stat, X.geyser, modele.geyser.continu.stat, M.step.geyser.continu.stat, trace.theta = TRUE, verbose = TRUE)


}





