# Fonction EM : spécifique à l'exemple choisi pour l'étape M
# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# it = iterations
EM <- function(theta, obs, modele.fun, M.step.fun, max.iter = 100, trace.theta = TRUE, epsilon = 1e-5){
  if(is.infinite(max.iter)) trace.theta <- FALSE
  l <- length(obs)

  if(trace.theta) {
    Theta <- matrix(NA_real_, ncol = max.iter, nrow = length(theta))
    rownames(Theta) <- names(theta)
    Theta[,1] <- theta
  }
 
  i <- 2
  repeat {
    mod <- modele.fun(theta, obs)
    if(any(is.infinite(mod$p.emiss))) {
      warning("Infinite density in model")
      i <- i-1
      break
    }
    # Etape E
    fo <- forward(mod)
    ba <- backward(fo)

    if(any(is.na(ba$phi))) {
      warning("Backward step failed")
      i <- i-1
      break
    }

    # Etape M
    theta1 <- M.step.fun(obs, ba)
    # when a latent state has probability 0, 
    # some of the parameters may be undefined: keep last value
    theta1 <- ifelse(is.na(theta1), theta, theta1)

    e <- sqrt( sum((theta1 - theta)**2) )
    theta <- theta1

    if(trace.theta) Theta[,i] <- theta
    if(e < epsilon | i == max.iter) break;
    i <- i+1
  }
  R <- list(theta = theta, iter = i)
  if(trace.theta) R$Theta <- Theta[, 1:i ]
  R
}


