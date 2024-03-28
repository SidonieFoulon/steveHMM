backward <- function(modele) {
  Tr <- modele$trans
  p.Em <- modele$p.emiss

  l <- ncol(p.Em)
  m <- nrow(Tr)

  alpha <- modele$alpha
  beta <- modele$beta
  delta <- array(NA_real_, dim = c(m, m, l), dimnames = list(rownames(Tr), rownames(Tr), NULL))

  phi <- matrix(NA_real_, nrow = m, ncol = l)
  rownames(phi) <- rownames(Tr)

  #initialisation du phi
  #delta_l n'a pas de valeur init
  phi[,l] <- beta[,l]

  #calcul du delta et phi simultane
  for(i in l:2){
    phia <- ifelse(alpha[,i] == 0, 0, phi[,i] / alpha[,i])
    delta[,,i] <- Tr * outer(beta[,i-1], phia, "*")
    phi[,i-1] <- rowSums(delta[,,i])
  }

  delta[,,1] <- rep(0,length(delta[,,1]))

  modele$delta <- delta
  modele$phi <- phi
  modele
}
