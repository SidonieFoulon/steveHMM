## on note : S in {"S1" : 1 ; "S2" : 2 ; "S3" : 3}
##           X in {"<=0" : 1 ; "1" : 2 ; "2" : 3}

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

  # proba emission
  f_emiss <- function(o) {
    pe <- rbind(ifelse(o == 1, 1-e21-e31, ifelse(o == 2, e21, e31)),
                ifelse(o == 1, e12, ifelse(o == 2, 1-e12-e32, e32)),
                ifelse(o == 1, e13, ifelse(o == 2, e23, 1-e13-e23)))
    rownames(pe) <- name.S
    pe
  }

  if(is.list(obs)) {
    p.emiss <- lapply(obs, f_emiss)
  } else {
    p.emiss <- f_emiss(obs)
  }

  # proba depart
  pi1 <- theta[13]
  pi2 <- theta[14]
  pi3 <- 1- pi1 - pi2
  pi <- c(s1 = pi1, s2 = pi2, s3 = pi3)

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.eye <- function(obs, backward) {

  ### transitions ###

  # first, expected number of times in each state (except last obs)
  f_state <- function(phi, i) {
    l <- ncol(phi)
    sum(phi[i, -l])
  }

  if( is.list(backward$phi) ) {
    p.S1 <- sum(sapply(backward$phi, f_state, 1))
    p.S2 <- sum(sapply(backward$phi, f_state, 2))
    p.S3 <- sum(sapply(backward$phi, f_state, 3))
  } else {
    p.S1 <- f_state(backward$phi, 1)
    p.S2 <- f_state(backward$phi, 2)
    p.S3 <- f_state(backward$phi, 3)
  }

  # then, expected number of transition from i to j
  f_trans <- function(delta) {
    rowSums(delta, dims = 2)
  }
  
  if( is.list(backward$delta) ) {
    D <- Reduce(`+`, lapply(backward$delta, f_trans))
  } else {
    D <- f_trans(backward$delta)
  }

  if(p.S1 > 0)
    t12 <- D[1,2] / p.S1
  else
    t12 <- 0

  if(p.S1 > 0)
    t13 <- D[1,3] / p.S1
  else
    t13 <- 0

  if(p.S2 > 0)
    t21 <- D[2,1] / p.S2
  else
    t21 <- 0

  if(p.S2 > 0)
    t23 <- D[2,3] / p.S2
  else
    t23 <- 0

  if(p.S3 > 0)
    t31 <- D[3,1] / p.S3
  else
    t31 <- 0

  if(p.S3 > 0)
    t32 <- D[3,2] / p.S3
  else
    t32 <- 0


  ### emissions ###
  # as before but including last state
  f_state <- function(phi, i) {
    sum(phi[i, ])
  }
  if( is.list(backward$phi) ) {
    p.S1 <- sum(sapply(backward$phi, f_state, 1))
    p.S2 <- sum(sapply(backward$phi, f_state, 2))
    p.S3 <- sum(sapply(backward$phi, f_state, 3))
  } else {
    p.S1 <- f_state(backward$phi, 1)
    p.S2 <- f_state(backward$phi, 2)
    p.S3 <- f_state(backward$phi, 3)
  }

  # emission of i in state s
  f_emiss <- function(phi, obs, i, s) {
    sum(phi[s, which(obs == i)])
  }
  if(is.list(backward$phi)) {
    e12 <- sum(mapply(f_emiss, backward$phi, obs, 1, 2)) / p.S2
    e13 <- sum(mapply(f_emiss, backward$phi, obs, 1, 3)) / p.S3
    e21 <- sum(mapply(f_emiss, backward$phi, obs, 2, 1)) / p.S1
    e23 <- sum(mapply(f_emiss, backward$phi, obs, 2, 3)) / p.S3
    e31 <- sum(mapply(f_emiss, backward$phi, obs, 3, 1)) / p.S1
    e32 <- sum(mapply(f_emiss, backward$phi, obs, 3, 2)) / p.S2
  } else {
    e12 <- f_emiss(backward$phi, obs, 1, 2) / p.S2
    e13 <- f_emiss(backward$phi, obs, 1, 3) / p.S3
    e21 <- f_emiss(backward$phi, obs, 2, 1) / p.S1
    e23 <- f_emiss(backward$phi, obs, 2, 3) / p.S3
    e31 <- f_emiss(backward$phi, obs, 3, 1) / p.S1
    e32 <- f_emiss(backward$phi, obs, 3, 2) / p.S2
  }

  ### etats initiaux ###
  if(is.list(backward$phi)) {
    pi1 <- mean(sapply(backward$phi, \(phi) phi[1,1]))
    pi2 <- mean(sapply(backward$phi, \(phi) phi[2,1]))
  } else {
    pi1 <- backward$phi[1,1]
    pi2 <- backward$phi[2,1]
  }

  c(t12, t13, t21, t23, t31, t32,
    e12, e13, e21, e23, e31, e32,
    pi1, pi2)
}

# version non "vectorisée"
M.step.eye.0 <- function(obs, backward) {
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

