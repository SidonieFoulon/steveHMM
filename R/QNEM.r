# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# it = iterations
QNEM <- function(theta, obs, modele.fun, M.step.fun, max.iter = 100, upper, lower, trace.theta = TRUE, reltol = sqrt(.Machine$double.eps), nb.em = length(theta), verbose = FALSE){
  c.armijo <- 1e-4
  tau_backt <- .5

  d <- length(theta)
  if(missing(lower)) lower <- rep(-Inf, d)
  if(missing(upper)) upper <- rep(-Inf, d)
  if(any(theta < lower | theta > upper)) stop("Initial value not in the bounding box")
  
  if(is.infinite(max.iter)) trace.theta <- FALSE
  if(trace.theta) {
    Theta <- matrix(NA_real_, ncol = max.iter, nrow = d)
    rownames(Theta) <- names(theta)
    Theta[,1] <- theta
  }
  i <- 2

  # keep trace of forward / backward calls
  nb.fw <- 0L
  nb.bw <- 0L

  # composantes bloquées
  I <- integer(0)
  # composantes non bloquées [on ne peut pas utiliser x[-I] quand I = integer(0) ...
  J <- seq_along(theta)

  # valeur initiale pour H
  # H <- diag(d)
  H <- matrix(NA_real_, d, d)

  # On amorce les calculs avant d'entrer dans la boucle
  # ** fo est supposé exister en début de boucle pour le point courant **
  # take modele derivatives
  mod <- modele_derivatives(modele.fun, theta, obs)
  fo <- forward_ll(mod, TRUE)
  nb.fw <- nb.fw + 1L
  ll <- fo$value
  gradient <- fo$gradient

  # the big loop
  repeat { 
    if(nb.em > 0) { # EM
      if(verbose) cat("EM step\n")
      # Finish the E step
      ba <- backward(mod, fo)
      nb.bw <- nb.bw + 1L

      # update H
      if(i > 2) {
        if(sum(s*y) > 0) H <- H_update(H, s, y)
      }
      # Etape M
      theta1 <- M.step.fun(obs, ba)

      # enforce bounds [rounding errors may occur event if the bounds are the natural ones]
      theta1 <- ifelse(theta1 > upper, upper, theta1)
      theta1 <- ifelse(theta1 < lower, lower, theta1)

      # in some cases, the M step produces NA's (e.g. conditionnal probabilities computed as 0/0)
      # this should be a correct way to deal with this
      theta1 <- ifelse(is.na(theta1), theta, theta1)

      # end of loop : update vars
      s <- theta1 - theta
      y <- fo$gradient - gradient
      theta <- theta1
      nb.em <- nb.em - 1
      # prépare la boucle suivante
      mod <- modele_derivatives(modele.fun, theta, obs)
      fo <- forward_ll(mod, TRUE)
      nb.fw <- nb.fw + 1L
      rel.ll <- abs(ll - fo$value) / (abs(ll) + reltol)
      ll <- fo$value
      gradient  <- fo$gradient
      QN <- FALSE
    } else { # QN
      if(verbose) cat("QN\n")
      if(any(is.na(H))) { # aucun pas d'EM n'a updaté H
        H <- diag(d)
      }
      # on peut arriver ici sans savoir que des variables sont bloqués
      # [ l'EM ne maintient pas I et J à jour ]
      over  <- (theta >= upper)
      under <- (theta <= lower)
      if(any(over[J] | under[J])) { # s'il y a des variables au bord non répertoriées
        # on commence par s'assurer qu'il n'y a pas de débordement
        theta <- ifelse(over, upper, theta)
        theta <- ifelse(under, lower, theta)
        # on met à jour la liste des variables bloquées 
        # [on regarde le gradient pour voir si c'est un vrai blocage]
        I1 <- which((over & gradient < 0) | (under & gradient > 0))
        I <- union(I, I1)
        J <- setdiff(J, I1)
        H <- restrict_inverse(H, I) # projection of H
      }
      # est-ce que des variables peuvent être débloqués ?
      I1 <- which( (over & gradient > 0) | (under & gradient) < 0)
      if(any(I1 %in% I)) {
        I <- setdiff(I, I1)
        J <- union(J, I1)
        H <- diag(d)
      }
      # if some variables are blocked, project the gradient
      gradient[I] <- 0
      if( sum(gradient**2) == 0 ) { # if it happens, the backtracking loop is infinitee
        break
      }

      # research direction
      p <- -(H %*% gradient)
      p.grad <- sum( (p * gradient)[J] )
      if(p.grad > 0) { # pas une bonne direction... ne devrait pas arriver
        stop("bad direction\n")
      }
      # preparing for backtracking
      lambda.i.max <- ifelse(p > 0, (upper - theta)/p, ifelse(p < 0, (lower - theta)/p, Inf))
      
      blocked.value <- ifelse(p >=0, upper, lower)
      lambda.max <- min(lambda.i.max[J], na.rm = TRUE)  # [J] : only unblocked vars
      if(lambda.max <= 0) { # ne devrait pas arriver
        stop("blocage non repéré")
      }
      lambda <- min(1, lambda.max)
      repeat { # bracktracking loop
        theta1 <- theta + lambda * p
        if(lambda == lambda.max) { # block variables
          # -> first take care of rounding errors 
          I1 <- which(lambda.i.max == lambda)
          theta1[I1] <- blocked.value[I1]
        }
        mod <- modele_derivatives(modele.fun, theta1, obs)
        fo1 <- forward_ll(mod, TRUE)
        nb.fw <- nb.fw + 1L
        ll1 <- if(is.na(fo1$value)) Inf else fo1$value
        gradient1 <- fo1$gradient
        gradient1[I] <- 0
        # Armijo rule
        if(ll1 <= ll + c.armijo * lambda * p.grad)
          break
        if(ll1 == ll) { # we may have converged...
          theta1 <- theta # don't move any more (will force an EM step)
          break
        }
        lambda <- lambda * tau_backt
        if(verbose) cat("backtracking\n")
      }
      # out of the backtracking loop with a "good" theta1
      # update H
      s <- theta1 - theta
      y <- gradient1 - gradient
      if(sum( (s*y)[J] ) <= 0) { # curvature condition [catches also the case theta1 = theta]
        if(verbose) cat("Non convex spot\n")
        nb.em <- 1
      } else {
        H <- H_update(H, s, y)
        # ------ if new wariables are blocked -----
        if(lambda == lambda.max) { 
          I <- union(I, I1)
          J <- setdiff(J, I1)
          if(verbose) cat("Blocking variables", I, "\n")
          gradient1[I] <- 0
          H <- restrict_inverse(H, I) # projection of H
        }
      }
      # update loop variables
      rel.ll <- abs(ll - ll1) / (abs(ll) + reltol)
      ll <- ll1
      fo <- fo1
      theta <- theta1
      gradient <- gradient1
      QN <- TRUE
    }
    if(trace.theta) Theta[,i] <- theta
    if(QN & nb.em > 0) { # in QN, decided to force at least one EM step 
      I <- integer(0)
      J <- seq_along(theta)
      H <- matrix(NA_real_, d, d)
    } else {
      if(rel.ll < reltol | i == max.iter) break 
    }
    i <- i+1
    if(verbose) cat("theta =", as.vector(theta), "\n")
    if(verbose) cat("neg log likelihood =", ll, "\n")
  }
  R <- list(theta = theta, iter = i, forwards = nb.fw, backwards = nb.bw)
  if(trace.theta) R$Theta <- Theta[, 1:i ]
  R
}

