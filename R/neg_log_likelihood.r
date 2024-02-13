# calcul de la vraissemblance à partir du forward de l'EM : probabilités conditionnelles
neg_log_likelihood <- function(theta, obs, modele.fun) {
  modele <- modele.fun(theta, obs)
  forward <- forward(modele)
  alpha <- forward$alpha
  p.Em <- modele$p.emiss
  -sum(log(colSums(alpha * p.Em)))
}


