forward <- function(modele) {
  Tr <- modele$trans
  p.Em <- modele$p.emiss

  l <- ncol(p.Em)
  m <- nrow(Tr)
 
  alpha <- matrix(NA_real_, ncol = l,nrow = m)
  rownames(alpha) <- rownames(Tr)
  
  beta <- matrix(NA_real_, ncol = l,nrow = m)
  rownames(beta) <- rownames(Tr)

  # initialisation
  # alpha1 = stationnaire
  alpha[,1] <- modele$pi

  tmp <- p.Em[,1]*alpha[,1]
  beta[,1] <- tmp/sum(tmp)
  
  for(i in 2:l){ 
    alpha[,i] <- beta[,i-1] %*% Tr
    
    tmp <- p.Em[,i]*alpha[,i]
    beta[,i] <- tmp/sum(tmp)
  }
  list(alpha = alpha, beta = beta) 
}

