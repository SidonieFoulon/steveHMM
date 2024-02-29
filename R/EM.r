# Fonction EM : spécifique à l'exemple choisi pour l'étape M
# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# it = iterations
EM <- function(theta, obs, modele.fun, M.step.fun, max.iter, epsilon = 1e-5){
  l <- length(obs)

  Theta <- matrix(NA_real_, ncol = max.iter, nrow = length(theta))
  rownames(Theta) <- names(theta)
  Theta[,1] <- theta

  for(i in 2:max.iter){
    mod <- modele.fun(theta, obs)

    # Etape E
    fo <- forward(mod)
    ba <- backward(mod, fo)

    # Etape M
    theta <- M.step.fun(obs, ba)

    Theta[,i] <- theta
    if( sqrt(sum((theta - Theta[,i-1])^2)) < epsilon)
      break;
  }
  Theta <- Theta[, apply(Theta, 2, function(x) !all(is.na(x)))]
  return(list(Theta = Theta, mod = mod, fb = c(fo, ba)))
}


