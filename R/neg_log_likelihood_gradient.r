neg_log_likelihood_gradient <- function(theta, obs, modele.fun) {
  modele <- modele.fun(theta, obs)

  Tr <- modele$trans
  p.Em <- modele$p.emiss

  d.p.Em <- modele$d.p.emiss
  d.Tr <- modele$d.trans

  l <- ncol(p.Em)
  m <- nrow(Tr)

  # nb de dérivées à prendre
  d <- length(modele$d.trans)

  # initialisation : alpha = état stationnaire
  alpha <- modele$pi
  d.alpha <- modele$d.pi
  d.tmp <- list()
  d.beta <- list()
  d.ll <- numeric(d)

  tmp <- p.Em[,1]*alpha
  for(k in 1:d) 
    d.tmp[[k]] <- d.p.Em[[k]][,1]*alpha + p.Em[,1]*d.alpha[[k]]

  stmp <- sum(tmp)
  beta <- tmp/stmp
  for(k in 1:d)
    d.beta[[k]] <- (stmp * d.tmp[[k]] - sum(d.tmp[[k]]) * tmp) / stmp**2

  # ll = log likelihood
  ll <- log(stmp)
  for(k in 1:d) 
    d.ll[k] <- sum(d.tmp[[k]]) / stmp
  
  for(i in 2:l){ 
    alpha <- beta %*% Tr
    for(k in 1:d)
      d.alpha[[k]] <- beta %*% d.Tr[[k]] + d.beta[[k]] %*%  Tr
     
    tmp <- p.Em[,i]*alpha
    for(k in 1:d) 
      d.tmp[[k]] <- d.p.Em[[k]][,i]*alpha + p.Em[,i]*d.alpha[[k]]

    stmp <- sum(tmp)
    beta <- tmp/stmp
    for(k in 1:d)
      d.beta[[k]] <- (stmp * d.tmp[[k]] - sum(d.tmp[[k]]) * tmp) / stmp**2

    ll <- ll + log(stmp)
    for(k in 1:d) 
      d.ll[k] <- d.ll[k] + sum(d.tmp[[k]]) / stmp
  }
  list(value = -ll, gradient = -d.ll)
}



