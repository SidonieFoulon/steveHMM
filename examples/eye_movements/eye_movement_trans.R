## on note : S in {"S1" : 1 ; "S2" : 2 ; "S3" : 3}
##           X in {"<=0" : 1 ; "1" : 2 ; "2" : 3}

p.stationnaire.eye <- function(pi1,pi2) {
  c( s1 = pi1, s2 = pi2, s3 = 1-pi1-pi2 )
}

modele.eye <- function(theta, obs, name.S = c("S1", "S2", "S3")) {
  # parametre de la matrice de transition
  t12 <- theta[1]
  t13 <- theta[2]
  t21 <- theta[3]
  t23 <- theta[4]
  t31 <- theta[5]
  t32 <- theta[6]


  trans <- matrix(c(1-t12-t13, t12, t13, t21, 1-t21-t23, t23, t31, t32, 1-t31-t32), nrow = 3, ncol = 3, byrow = TRUE)
  colnames(trans) <- rownames(trans) <- name.S

  if(any(trans < -1e-10)) trans[] <- NA_real_

  # parametre de l'émission
  e12 <- theta[7]
  e13 <- theta[8]
  e21 <- theta[9]
  e23 <- theta[10]
  e31 <- theta[11]
  e32 <- theta[12]

  p.emiss <- rbind(ifelse(obs == 1, 1-e21-e31, ifelse(obs == 2, e21, e31)),
                   ifelse(obs == 1, e12, ifelse(obs == 2, 1-e12-e32, e32)),
                   ifelse(obs == 1, e13, ifelse(obs == 2, e23, 1-e13-e23)))
  rownames(p.emiss) <- name.S

  # proba
  pi1 <- theta[13]
  pi2 <- theta[14]
  pi3 <- 1- pi1 - pi2
  pi <- c(s1 = pi1, s2 = pi2, s3 = pi3)

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.eye <- function(obs, backward) {
  l <- ncol(backward$phi)


  ### transitions ###
  p.S1 <- sum(backward$phi[1,-l])
  p.S2 <- sum(backward$phi[2,-l])
  p.S3 <- sum(backward$phi[3,-l])

  D12 <- sum(backward$delta[1,2,-1])
  if(p.S1 > 0)
    t12 <- D12 / p.S1
  else
    t12 <- 0

  D13 <- sum(backward$delta[1,3,-1])
  if(p.S1 > 0)
    t13 <- D13 / p.S1
  else
    t13 <- 0

  D21 <- sum(backward$delta[2,1,-1])
  if(p.S2 > 0)
    t21 <- D21 / p.S2
  else
    t21 <- 0

  D23 <- sum(backward$delta[2,3,-1])
  if(p.S2 > 0)
    t23 <- D23 / p.S2
  else
    t23 <- 0

  D31 <- sum(backward$delta[3,1,-1])
  if(p.S3 > 0)
    t31 <- D31 / p.S3
  else
    t31 <- 0

  D32 <- sum(backward$delta[3,2,-1])
  if(p.S3 > 0)
    t32 <- D32 / p.S3
  else
    t32 <- 0


  ### emissions ###
  e12 <- sum(backward$phi[2,which(obs==1)]) / sum(backward$phi[2,])
  e13 <- sum(backward$phi[3,which(obs==1)]) / sum(backward$phi[3,])
  e21 <- sum(backward$phi[1,which(obs==2)]) / sum(backward$phi[1,])
  e23 <- sum(backward$phi[3,which(obs==2)]) / sum(backward$phi[3,])
  e31 <- sum(backward$phi[1,which(obs==3)]) / sum(backward$phi[1,])
  e32 <- sum(backward$phi[2,which(obs==3)]) / sum(backward$phi[2,])


  ### etats initiaux ###
  pi1 <- backward$phi[1,1]
  pi2 <- backward$phi[2,1]

  c(t12, t13, t21, t23, t31, t32,
    e12, e13, e21, e23, e31, e32,
    pi1, pi2)
}

# nos observations  :
X.eye <- read.csv("/home/sidonie/Téléchargements/em-y35-fasttext.csv")
X.eye <- X.eye[,c(2,4,18)]
X.eye$READMODE[which(X.eye$READMODE < 0)] <- 0
X.eye$READMODE <- X.eye$READMODE +1

# nos paramètres d'initialisation :
pi.eye <- p.stationnaire.eye(pi1 = 0.1, pi2 = 0.2)
par.eye <- c(  t12 <- 0.2,
               t13 <- 0.1,
               t21 <- 0.2,
               t23 <- 0.2,
               t31 <- 0.1,
               t32 <- 0.2,
               e12 <- 0.2,
               e13 <- 0.1,
               e21 <- 0.2,
               e23 <- 0.2,
               e31 <- 0.1,
               e32 <- 0.2,
               pi.eye[1],
               pi.eye[2])


if(FALSE) {
  library(steveHMM)

  # qN
  qn.dich <- quasi_newton(par.dich, X.dich, modele.geyser, lower = rep(0,7) + 0.01, upper = rep(1,7) - 0.01)
  qn.dich.trace <- capture_quasi_newton(par.dich, X.dich, modele.geyser, lower = rep(0,7) + 0.01, upper = rep(1,7) - 0.01)

  # EM
  em.eye <- EM(par.eye, X.eye$READMODE[1:29], modele.eye, M.step.eye, max.iter = 100, trace.theta = TRUE, criteria = "rel")

  # SQUAREM
  squarem.eye <- SQUAREM(par.eye, X.eye$READMODE[1:29], modele.eye, M.step.eye, lower = rep(0,14), upper = rep(1,14), max.iter = 1000, trace.theta = TRUE, criteria = "reltol")

  # QNEM
  qnem.eye <- QNEM(par.eye, X.eye$READMODE[1:29], modele.eye, M.step.eye, max.iter = 1000, lower = rep(0,14), upper = rep(1,14), trace.theta = TRUE, verbose = TRUE)
}
