# Fonction EM : spécifique à l'exemple choisi pour l'étape M
# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# it = iterations
EM <- function(theta, obs, modele.fun, M.step.fun, it){
  l <- length(obs)
  
  Theta <- matrix(NA_real_, ncol = it, nrow = length(theta))
  rownames(Theta) <- names(theta)
  Theta[,1] <- theta

  for(i in 2:it){
    mod <- modele.fun(theta, obs) 

    # Etape E
    fo <- forward(mod)
    ba <- backward(mod, fo)
  
    # Etape M 
    theta <- M.step.fun(obs, ba)

    Theta[,i] <- theta
  }
  return(list(Theta = Theta, mod = mod, fb = c(fo, ba)))
}


