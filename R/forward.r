#' Forward algorithm
#'
#' @param modele a model function
#'
#' @return A list with components `alpha`, `beta`, and `likelihood`. `alpha` and `beta` are matrices with the
#' the forward probabilities, a and `likelihood` is the log-likelihood.
#'
#' @export
#'


forward <- function(modele) {
  Tr <- modele$trans
  p.Em <- modele$p.emiss
  pi <- modele$pi

  if(is.list(p.Em)) {
    n <- length(p.Em)
    mo <- modele

    modele$alpha <- vector("list", n)
    modele$beta  <- vector("list", n)
    ll <- 0
    for(i in 1:n) {
      mo$p.emiss <- p.Em[[i]]
      fo <- forward(mo)
      modele$alpha[[i]] <- fo$alpha
      modele$beta[[i]]  <- fo$beta
      ll <- ll + fo$likelihood
    }
    modele$likelihood <- ll
    return(modele)
  }

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

