library(steveHMM)
library(salad)

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

test.time.gc <- function(theta, obs, it){

  #qN
  cpt <- 0
  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    X <- sample(obs)
    debut <- Sys.time()
    try(quasi_newton(theta, X, modele.geyser.continu, trace = FALSE, lower = rep(0.01,8), upper = c(0.99,0.99, rep(Inf,6))))
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
    try(EM(theta, X, modele.geyser.continu, M.step.geyser.continu, max.iter = Inf, trace.theta = FALSE))
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
    try(SQUAREM(theta, X, modele.geyser.continu, M.step.geyser.continu, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = FALSE))
    fin <- Sys.time()
    cpt <- cpt + (fin-debut)
  }

  SQ <- cpt



  return(list(BW = BW , SQUAREM = SQ, qN = qN))
}

#observations
library(MASS)
X.geyser <- faithful$eruptions
par.geyser <- c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39)

#test du temps
set.seed(28)
tps_gc <- test.time.gc(par.geyser, obs = X.geyser, it = 5000)
tps_gc
saveRDS(tps_gc, "/home/sidonie/Bureau/github/steveHMM/examples/tps_geyser_continu.rds")
