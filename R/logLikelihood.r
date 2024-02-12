logLikelihood <- function(theta, obs, modele.fun) {
  modele <- modele.fun(theta, obs)
  forward <- forward(modele)
  alpha <- forward$alpha 
  p.Em <- modele$p.emiss
  sum(log(colSums(alpha * p.Em)))   
}
