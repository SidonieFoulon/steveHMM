# fonction qui prend un theta numérique et appelle 'dual' seulement
# pour obtenir les dérivées des matrices envoyées par modèle
# ensuite la differenciation automatique est faite par la fonction elle même
neg_log_likelihood_gradient <- function(theta, obs, modele.fun) {

  theta1 <- dual(theta)  
  vn <- varnames(theta1)

  modele <- modele.fun(theta1, obs)

  Tr <- value(modele$trans)
  p.Em <- value(modele$p.emiss)
  Pi <- value(modele$pi)

  d.Tr <- d(modele$trans, vn)
  d.p.Em <- d(modele$p.emiss, vn)
  d.Pi <- d(modele$pi, vn)

  l <- ncol(p.Em)
  m <- nrow(Tr)

  # nb de dérivées à prendre
  d <- length(vn)

  # initialisation : alpha = état stationnaire
  alpha <- Pi
  d.alpha <- d.Pi
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



