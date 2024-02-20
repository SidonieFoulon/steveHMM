# calcul de la log vraisemblance à partir du forward 'classique'
# On récrit l'algo forward : cette version ne conserve pas les valeurs calculées
# (on ne fait pas un backward ensuite)
neg_log_likelihood <- function(theta, obs, modele.fun) {
  modele <- modele.fun(theta, obs)
  Tr <- modele$trans
  p.Em <- modele$p.emiss
  l <- ncol(p.Em)

  # initialisation : alpha = état stationnaire
  alpha <- modele$pi
  tmp <- p.Em[,1]*alpha
  stmp <- sum(tmp)
  beta <- tmp/stmp

  # ll = log likelihood
  ll <- log(stmp)

  for(i in 2:l){
    alpha <- beta %*% Tr
    tmp <- p.Em[,i]*alpha
    stmp <- sum(tmp)
    beta <- tmp/stmp
    ll <- ll + log(stmp)
  }
  -ll
}


if(FALSE) { # vieille version, à partir de la fonction forward 
neg_log_likelihood <- function(theta, obs, modele.fun) {
  modele <- modele.fun(theta, obs)
  forward <- forward(modele)
  alpha <- forward$alpha
  p.Em <- modele$p.emiss
  -sum(log(colSums(alpha * p.Em)))
}
}
