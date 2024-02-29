# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# upper, lower = comme dans optim
# max.iter = max iterations
# epsilon = critere de convergence
SQUAREM <- function(theta, obs, modele.fun, M.step.fun, lower, upper, max.iter, epsilon = 1e-5){
  l <- length(obs)
  if(missing(lower)) lower <- rep(-Inf, length(theta))
  if(missing(upper)) upper <- rep(+Inf, length(theta))
  if( any(theta < lower) | any(theta > upper) )
    stop("Bad starting point")

  # to keep all iterates
  Theta <- matrix(NA_real_, ncol = max.iter, nrow = length(theta))
  rownames(Theta) <- names(theta)
  Theta[,1] <- theta

  # U is for stacking a point + two EM itererates
  U <- matrix(NA_real_, ncol = 3, nrow = length(theta))
  U[,1] <- theta

  # 1st EM iterate need to be computed before entering the loop
  mod <- modele.fun(theta, obs)
  fo <- forward(mod)
  ba <- backward(mod, fo)
  theta <- M.step.fun(obs, ba)
  U[,2] <- theta

  Theta[,2] <- theta
  k <- 3


  repeat { # The big loop

    repeat { # iterate EM until beta > 0
      mod <- modele.fun(theta, obs)
      fo <- forward(mod)
      ba <- backward(mod, fo)
      theta <- M.step.fun(obs, ba)
      U[,3] <- theta
      # *** keep track of theta ***
      Theta[,k] <- theta
      k <- k+1;
      if(k > max.iter){
        Theta <- Theta[, apply(Theta, 2, function(x) !all(is.na(x)))]
        return(Theta)
      }

      # check convergence
      if( sqrt(sum((U[,2] - U[,3])**2)) < epsilon ){
        Theta <- Theta[, apply(Theta, 2, function(x) !all(is.na(x)))]
        return(Theta)
      }


      beta <- squarem.beta(U)
      if(beta <= 0) { # continue EM iterations
        U[,1:2] <- U[,2:3] # shift U
      } else break; # we got a beta > 0
    }

    # computing log likelihood
    # first for the last theta
    ll0 <- sum(log(colSums(fo$alpha * mod$p.emiss)))

    # capping beta
    if(beta > 2) beta <- 2
    # now for the potential new theta, until it's higher than ll0
    repeat { # backtracking loop
      theta1 <- squarem.proposal(U, beta)
      # checking box bounds
      if( any(theta1 < lower) | any(theta1 > upper) ) {
        # a good idea could be to go on the box boundary
        # but for now backtracking is ok
        ll1 <- -Inf
      } else {
        mod <- modele.fun(theta1, obs)
        fo <- forward(mod)
        ll1 <- sum(log(colSums(fo$alpha * mod$p.emiss)))
      }
      if(ll1 >= ll0) { # accept proposal
        theta <- theta1
        U[,1] <- theta # it's our new starting point
        # take advantage that the forward has been done already to finish the E step
        ba <- backward(mod, fo)
        # *** keep track of theta ***
        Theta[,k] <- theta
        k <- k+1;
        if(k > max.iter){
          Theta <- Theta[, apply(Theta, 2, function(x) !all(is.na(x)))]
          return(Theta)
          }

        # M step
        theta <- M.step.fun(obs, ba)
        U[,2] <- theta
        # *** keep track of theta ***
        Theta[,k] <- theta
        k <- k+1;
        if(k > max.iter){
          Theta <- Theta[, apply(Theta, 2, function(x) !all(is.na(x)))]
          return(Theta)
        }

        # check convergence
        if( sqrt(sum((U[,1] - U[,2])**2)) < epsilon ){
          Theta <- Theta[, apply(Theta, 2, function(x) !all(is.na(x)))]
          return(Theta)
        }



        break # break from the backtracking loop, ready to compute U[,3]
      }
      # else, backtracking
      beta <- beta/2
# cat("backtracking beta =", beta, "\n")
      if(beta <= 0.25) { # we bactracked too far, going back to EM
        beta <- 0
        U[,1:2] <- U[,2:3] # shift U, readdy for going back to beginning of "big loop"
        break
      }
    }
  }

}


squarem.beta <- function(U) {
  r <- U[,2] - U[,1]
  v <- U[,3] - 2*U[,2] +  U[,1]
  beta <- sqrt(sum(r**2)) / sqrt(sum(v**2)) - 1
  # cat("beta = ", beta, "\n")
  beta
}

squarem.proposal <- function(U, beta) {
  if(beta <= 0) return(U[,3])
  U[,2] + 2*beta*(U[,3] - U[,2]) + beta^2*(U[,3] - 2*U[,2] +  U[,1])
}
