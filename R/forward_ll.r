# forward that compute neg log likelihood and its gradient
forward_ll <- function(modele, keep.forward) {

  Tr <- modele$trans
  p.Em <- modele$p.emiss
  Pi <- modele$pi

  d.Tr <- modele$dtrans
  d.p.Em <- modele$dp.emiss
  d.Pi <- modele$dpi

  if(is.list(p.Em)) {
    n <- length(p.Em)
    mo <- modele

    if(keep.forward) {
      modele$alpha <- vector("list", n)
      modele$beta  <- vector("list", n)
    }

    ll <- 0
    gr <- 0
    for(i in 1:n) {
      mo$p.emiss <- p.Em[[i]]
      mo$dp.emiss <- d.p.Em[[i]]
      fo <- forward_ll(mo, keep.forward)
      if(keep.forward) {
        modele$alpha[[i]] <- fo$alpha
        modele$beta[[i]]  <- fo$beta
      }
      ll <- ll + fo$likelihood
      gr <- gr + fo$likelihood.gradient
    }
    modele$likelihood <- ll
    modele$likelihood.gradient <- gr
    return(modele)
  }


  l <- ncol(p.Em)
  m <- nrow(Tr)

  # matrices pour garder alpha, beta, si keep.forward est TRUE
  if(keep.forward) {
    ALPHA <- matrix(NA_real_, ncol = l,nrow = m)
    rownames(ALPHA) <- rownames(Tr)

    BETA <- matrix(NA_real_, ncol = l,nrow = m)
    rownames(BETA) <- rownames(Tr)
  }

  # nb de dérivées à prendre
  d <- length(d.Tr)

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
  
  # trace alpha, bata
  if(keep.forward) {
    ALPHA[,1] <- alpha
    BETA[,1] <- beta
  }

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

    # trace alpha, beta
    if(keep.forward) {
      ALPHA[,i] <- alpha
      BETA[,i] <- beta
    }
  }

  modele$likelihood <- -ll
  modele$likelihood.gradient <- -d.ll
  if(keep.forward) {
    modele$alpha <- ALPHA
    modele$beta <- BETA
  }
  modele
}

