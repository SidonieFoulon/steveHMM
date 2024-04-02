BFGS <- function(theta, FUN, max.iter = 100, upper, lower, trace.theta = TRUE, epsilon = 1e-5){
  nb.em <- 10
  c.armijo <- 1e-4
  tau_backt <- .5

  H <- diag(length(theta))

  # the big loop
  repeat {
      fo <- FUN(theta)
      ll <- fo$likelihood
      gradient <- fo$likelihood.gradient
      # research direction
      p <- -(H %*% gradient)
      p.grad <- sum(p * gradient)
      if(p.grad > 0) { # pas une bonne direction...
        stop("bad direction")
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
        fo1 <- FUN(theta1)
        ll1 <- if(is.na(fo1$likelihood)) Inf else fo1$likelihood
        # Armijo rule [IL FAUT AUSSI VERIFIER LA CURVATURE RULE]
        if(ll1 < ll + c.armijo * lambda * p.grad)
          break
        lambda <- lambda * tau_backt
      }
      # out of the backtracking loop with a "good" theta1
      # update H
      s <- theta1 - theta
      y <- fo1$likelihood.gradient - gradient
      if(sum(s*y) <= 0) { # curvature condition 
        stop("curvature condition does not hold*\n")
      }
      H <- H_update(H, s, y)
      theta <- theta1
      gradient <- fo1$likelihood.gradient 
    print(H)
    print(theta)
  }
}

