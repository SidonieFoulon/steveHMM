# theta = (a, b, e)
# a et b = transition probas, e = proba emission error
# obs = vecteur de 1 et 2
modele.2s <- function(theta, obs) {
  a <- theta[1]
  b <- theta[2]
  trans <- matrix( c(1-a, a, b, 1-b), byrow = TRUE, nrow = 2, ncol = 2)
  colnames(trans) <- rownames(trans) <- c("S1", "S2")
  # etat stationnaire
  stat <- c(b, a)/(a + b)
  # proba d'erreur d'émission
  e <- theta[3]
  p.emiss <- rbind(S1 = ifelse(obs == 1, 1-e, e), S2 = ifelse(obs == 1, e, 1-e))
  list(trans = trans, pi = stat, p.emiss = p.emiss, theta = theta)
}

M.step.2s <- function(obs, backward) {
  D11 <- sum(backward$delta[1,1,])
  D12 <- sum(backward$delta[1,2,]) + backward$phi[2,1]
  D21 <- sum(backward$delta[2,1,]) + backward$phi[1,1]
  D22 <- sum(backward$delta[2,2,])
  a <- backward$theta[1]
  b <- backward$theta[2]
  a <- (D12 - a*(1-a)/(a+b)) / (D11 + D12)
  b <- (D21 - b*(1-b)/(a+b)) / (D21 + D22)

  e <- (sum(backward$phi[2,which(obs == 1)]) + sum(backward$phi[1,which(obs == 2)])) / ncol(backward$phi)
  c(a, b, e)
}

X.2s <- rep(rep(1:2, 115), c(3,1,6,1,2,1,3,3,13,1,15,1,3,1,7,
1,3,1,2,1,5,1,4,9,1,10,1,3,5,1,1,10,5,1,12,16,1,9,1,1,1,12,5,
1,6,1,3,1,6,1,7,3,2,3,10,1,15,3,6,1,4,1,8,1,11,1,12,7,4,15,5,
2,5,7,3,14,5,1,1,1,14,1,6,6,4,1,9,1,6,6,1,1,1,1,1,1,1,4,7,1,
2,1,6,1,5,7,14,12,1,1,11,1,1,1,6,7,1,10,1,1,3,3,3,1,6,6,10,1,
11,6,1,7,5,13,2,2,3,2,17,1,2,1,18,1,4,1,7,1,2,1,6,1,7,1,3,1,
2,1,1,1,8,2,15,1,7,1,3,1,2,1,9,2,16,12,1,2,1,1,4,1,4,1,4,1,
12,6,5,3,1,1,10,1,2,1,7,1,2,2,9,4,1,12,9,2,1,10,1,1,1,22,1,8,
1,4,1,1,2,1,11,1,8,1,11,1,2,1,1,9,6,12))

par.2s <- c(0.1, 0.1, 0.2)


if(FALSE) {
require(steveHMM)
qnem <- QNEM(par.2s, X.2s, modele.2s, M.step.2s, lower = c(0,0,0), upper = c(1,1,1), verbose = TRUE)
em <- EM(par.2s, X.2s, modele.2s, M.step.2s)
}
