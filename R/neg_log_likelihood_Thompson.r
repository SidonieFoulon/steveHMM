# autre proposition pour le calcul de la vraisemblance : probabilités jointes (Thompson ?)

neg_log_likelihood_Thompson <- function(theta, obs, modele.fun) {

  # Appel du modele
  modele <- modele.fun(theta, obs)

  # Calcul de la vraisemblance par un forward
  L <- length(obs)
  alpha <- matrix(0, ncol = L, nrow = nrow(modele$trans))
  rownames(alpha) <- rownames(modele$trans)

  # Initialisation
  alpha[,1] <- modele$pi

  # Recurrence
  for(i in 2:L){
    alpha[,i] <- ( modele$p.emiss[,i-1] * alpha[,i-1] )%*% modele$trans
  }

  # Vraisemblance
  likelihood <- sum(modele$p.emiss[,L] * alpha[, L]) #la dernière colonne de la chaine contient toute l'information de la chaine

  # - log vraisemblance (a minimiser)
  return(-log(likelihood))
}
