forward <- function(modele) {
  Tr <- modele$trans
  p.Em <- modele$p.emiss
  pi <- modele$pi

  dd <- is(Tr, "dual") | is(p.Em, "dual") | is(pi, "dual")
  if(dd) vn <- unique(c(varnames(Tr), varnames(p.Em), varnames(pi)))

  l <- ncol(p.Em)
  m <- nrow(Tr)

  alpha <- matrix(NA_real_, ncol = l,nrow = m)
  rownames(alpha) <- rownames(Tr)

  beta <- matrix(NA_real_, ncol = l,nrow = m)
  rownames(beta) <- rownames(Tr)

  if(dd) {
    alpha <- dual(alpha, varnames = vn, constant = TRUE)
    beta  <- dual(beta,  varnames = vn, constant = TRUE)
  }

  # initialisation
  # alpha1 = stationnaire
  alpha[,1] <- pi

  tmp <- p.Em[,1]*alpha[,1]
  stmp <- sum(tmp)
  beta[,1] <- tmp/stmp
  ll <- log(stmp)

  for(i in 2:l){
    alpha[,i] <- beta[,i-1] %*% Tr

    tmp <- p.Em[,i]*alpha[,i]
    stmp <- sum(tmp)
    beta[,i] <- tmp/stmp
    ll <- ll + log(stmp)
  }
  modele$alpha <- alpha
  modele$beta <- beta
  modele$likelihood <- ll
  modele
}

