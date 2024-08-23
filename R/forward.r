<<<<<<< HEAD
#' Forward algorithm
#'
#' @param modele a model function
#'
#' @export
#'
=======
>>>>>>> fbedfc2faf19b83b437296c93248a479e48507e4

forward <- function(modele) {
  Tr <- modele$trans
  p.Em <- modele$p.emiss
  pi <- modele$pi

  l <- ncol(p.Em)
  m <- nrow(Tr)

  alpha <- matrix(NA_real_, ncol = l,nrow = m)
  rownames(alpha) <- rownames(Tr)

  beta <- matrix(NA_real_, ncol = l,nrow = m)
  rownames(beta) <- rownames(Tr)

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

