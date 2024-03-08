# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# it = iterations
QNEM <- function(theta, obs, modele.fun, M.step.fun, max.iter = 100, upper, lower, trace.theta = TRUE, epsilon = 1e-5){
  if(is.infinite(max.iter)) trace.theta <- FALSE
  l <- length(obs)

  if(trace.theta) {
    Theta <- matrix(NA_real_, ncol = max.iter, nrow = length(theta))
    rownames(Theta) <- names(theta)
    Theta[,1] <- theta
  }
 
  i <- 2
  nb.em <- 10
  c.armijo <- 1e-4
  tau_backt <- .5

  H <- diag(length(theta))

     
  # On amorce les calculs avant d'entrer dans la boucle
  # ** fo est supposé exister en début de boucle pour le point courant **
  # take modele derivatives
  mod <- modele_derivatives(modele.fun, theta, obs)
  # étape E
  fo <- forward_ll(mod, TRUE)
  # the big loop
  repeat {
    if(nb.em > 0) { # EM
cat("EM\n")
      # take modele derivatives
      # mod <- modele_derivatives(modele.fun, theta, obs)

      # fin étape E
      # fo <- forward_ll(mod, TRUE)
      ba <- backward(mod, fo)
      
      ll <- fo$value
      # update H
      # s = theta[i-1] - theta[i-2]
      # gradient = gradient[i - 2]
      # fo$gradient = gradient[i - 1]
       if(i > 2) {
        y <- fo$gradient - gradient
        if(sum(s*y) > 0) H <- H_update(H, s, y)
      }

      # Etape M
      theta1 <- M.step.fun(obs, ba)

      # update vars
      e <- sqrt( sum((theta1 - theta)**2) )
      s <- theta1 - theta
      theta <- theta1
      gradient <- fo$gradient
      nb.em <- nb.em - 1
      # prépare la boucle suivante
      mod <- modele_derivatives(modele.fun, theta, obs)
      fo <- forward_ll(mod, TRUE)
    } else { # QN 
browser()
cat("QN\n")
      # mod <- modele_derivatives(modele.fun, theta, obs)
      # compute gradient at current point
      # fo <- forward_ll(mod, TRUE)
      ll <- fo$value
      gradient <- fo$gradient
      # research direction
      p <- -(H %*% gradient)
      p.grad <- sum(p * gradient)
      if(p.grad > 0) { # pas une bonne direction...
stop("bad direction\n")
        nb.em <- 1
        next
      }
      # preparing for backtracking
      lambda.i.max <- ifelse(p >= 0, (upper - theta)/p, (lower - theta)/p)
      blocked.value <- ifelse(p >=0, upper, lower)
      lambda.max <- min(lambda.i.max)
      lambda <- min(1, lambda.i.max)
      repeat { # bracktracking loop
        theta1 <- theta + lambda * p
        if(lambda == lambda.max) { # take care of rounding errors
          I <- which(lambda.i.max == lambda.max)
          theta1[ I ] <- blocked.value[I]
        } 
        mod <- modele_derivatives(modele.fun, theta1, obs)
        fo1 <- forward_ll(mod, TRUE)
        ll1 <- if(is.na(fo1$value)) Inf else fo1$value
        # Armijo rule
        if(ll1 < ll + c.armijo * lambda * p.grad)
          break
        lambda <- lambda * tau_backt
      }
      # out of the backtracking loop with a "good" theta1
      # update H
      s <- theta1 - theta
      y <- fo1$gradient - gradient
      if(sum(s*y) <= 0) { # curvature condition 
        # GO BACK TO EM ?
        browser();
        fo <- fo1
        theta <- theta1
        gradient <- fo1$gradient
        nb.em <- 1
        next
      }
      # update H, theta, gradient, etc
      H <- H_update(H, s, y)
      fo <- fo1
      theta <- theta1
      gradient <- fo1$gradient
      e <- sqrt( sum(s**2) )
    }
    if(trace.theta) Theta[,i] <- theta
    if(e < epsilon | i == max.iter) break;
    i <- i+1
#print(H)
#print(theta)
cat("neg log likelihood =", ll, "\n")
  }
  R <- list(theta = theta, iter = i)
  if(trace.theta) R$Theta <- Theta[, 1:i ]
  R
}

