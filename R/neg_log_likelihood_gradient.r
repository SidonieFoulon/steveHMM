# fonction qui prend un theta numérique et appelle 'dual' seulement
# pour obtenir les dérivées des matrices envoyées par modèle
# ensuite la differenciation automatique est faite par la fonction elle même
neg_log_likelihood_gradient <- function(theta, obs, modele.fun) {
  mo <- modele_derivatives(modele.fun, theta, obs, modele.fun)
  forward_ll(mo, keep.forward = FALSE)
}

