library(steveHMM)
library(salad)

path <- "/home/sidonie/Bureau/github/steveHMM/examples/"

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
    a <- (D13 - a0*(1-a0)*b0/(2 - 2*b0 + a0*b0)) / (D12 + D13)
  else
    a <- 1 # pas d'état court

  if(D31 + D33 > 0)
    b <- (D33 - b0*(1-b0)*(a0-2)/(2 - 2*b0 + a0*b0)) / (D31 + D33)
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

  r <- c(a = as.numeric(a), b = as.numeric(b), muc = muc, mul = mul, muls = muls, sdc = sdc, sdl = sdl, sdls = sdls)
  r
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

  # etat stationnaire
  pi_c <- theta[9]
  pi_l <- theta[10]
  pi_ls <- 1 - pi_c - pi_l
  pi <- c( c = pi_c, l = pi_l, ls = pi_ls )

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.geyser.continu <- function(obs, backward) {
  l <- ncol(backward$phi)

  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  D13 <- sum(backward$delta[1,3,-1])
  p.court <- sum(backward$phi[1,-l]) # = D12 + D13
  if(p.court > 0)
    a <- D13 / p.court
  else
    a <- 1 # pas d'état court

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  D33 <- sum(backward$delta[3,3,-1])
  p.long.st <- sum(backward$phi[3,-l]) # = D31 + D33
  if(p.long.st > 0)
    b <- D33 / p.long.st
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

  #etats initiaux
  pi_c <- backward$phi[1,1]
  pi_l <- backward$phi[2,1]

  r <- c(a = a, b = b, muc = muc, mul = mul, muls = muls, sdc = sdc, sdl = sdl, sdls = sdls, pi_c = pi_c, pi_l = pi_l)
  r
}

test.init.gc <- function(obs, it){
  qN <- matrix(nrow = it, ncol = 22)
  colnames(qN) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l", "tps")
  BW <- matrix(nrow = it, ncol = 22)
  colnames(BW) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l", "tps")
  SQ <- matrix(nrow = it, ncol = 24)
  colnames(SQ) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l","forward", "backward", "tps")
  QNEM <- matrix(nrow = it, ncol = 24)
  colnames(QNEM) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l", "forward", "backward", "tps")

  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    par <- c(runif(2), runif(3, 0, 6), runif(3,1,3))

    d <- Sys.time()
    tr <- try(qn <- quasi_newton(par, obs, modele.geyser.continu.stat, lower = rep(0.01,8), upper = c(0.99,0.99, rep(Inf,6))))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      qN[i,1] <- qn$counts[1]
      qN[i,2:9] <- qn$par
      qN[i,12:19] <- par
      qN[i,22] <- f-d
    }

    test.em <- function(par){
      em <- EM(par, obs, modele.geyser.continu.stat, M.step.geyser.continu.stat, max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(em)}
    d <- Sys.time()
    tr <- try(em <- test.em(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      BW[i,1] <- em$iter
      BW[i,2:9] <- em$theta
      BW[i,12:19] <- par
      BW[i,22] <- f-d
    }

    test.sq <- function(par){
      squarem <- SQUAREM(par, obs, modele.geyser.continu.stat, M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(squarem)}
    d <- Sys.time()
    tr <- try(squarem <- test.sq(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      SQ[i,1] <- squarem$iter
      SQ[i,2:9] <- squarem$theta
      SQ[i,12:19] <- par
      SQ[i,22] <- squarem$forwards
      SQ[i,23] <- squarem$backwards
      SQ[i,24] <- f-d
    }

    test.qnem <- function(par){
      qnem <- QNEM(par, obs, modele.geyser.continu.stat, M.step.geyser.continu.stat )
      return(qnem)}
    d <- Sys.time()
    tr <- try(qnem <- test.qnem(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      QNEM[i,1] <- qnem$iter
      QNEM[i,2:9] <- qnem$theta
      QNEM[i,12:19] <- par
      QNEM[i,22] <- qnem$forwards
      QNEM[i,23] <- qnem$backwards
      QNEM[i,24] <- f-d
    }
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN, QNEM = QNEM))
}


#observations
library(MASS)
X.geyser <- faithful$eruptions

#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it_gc <- test.init.gc(X.geyser, 1000)
saveRDS(nb_it_gc, paste0(path,"nb_it_geyser_continu_qnem_meot.rds"))
