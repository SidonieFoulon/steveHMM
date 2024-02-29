# theta = paramètres de début
# obs = observations
# modele.fun = une fonction "modèle"
# M.step.fun = une fonction M step (input : obs, backward, output : theta)
# upper, lower = comme dans optim
# max.iter = max iterations
# epsilon = critere de convergence
SQUAREM <- function(theta, obs, modele.fun, M.step.fun, lower, upper, max.iter = 100, trace.theta = TRUE, epsilon = 1e-5) {
  if(is.infinite(max.iter)) trace.theta <- FALSE

  l <- length(obs)
  if(missing(lower)) lower <- rep(-Inf, length(theta))
  if(missing(upper)) upper <- rep(+Inf, length(theta))
  if( any(theta < lower) | any(theta > upper) ) 
    stop("Bad starting point")

  # to keep all iterates
  if(trace.theta) {
    Theta <- matrix(NA_real_, ncol = max.iter, nrow = length(theta))
    rownames(Theta) <- names(theta)
    Theta[,1] <- theta
  }

  # keep trace of forward / backward calls
  nb.fw <- 0L
  nb.bw <- 0L

  # U is for stacking a point + two EM itererates
  U <- matrix(NA_real_, ncol = 3, nrow = length(theta))
  U[,1] <- theta

  # 1st EM iterate need to be computed before entering the loop
  mod <- modele.fun(theta, obs) 
  fo <- forward(mod)      
  ba <- backward(mod, fo)  
  nb.fw <- nb.fw + 1
  nb.bw <- nb.bw + 1

  theta <- M.step.fun(obs, ba)
  U[,2] <- theta
  # *** keep trace of theta ***
  if(trace.theta) Theta[,2] <- theta
  k <- 3
 
  EXIT <- FALSE 
  repeat { # The big loop

    repeat { # iterate EM until beta > 0
      mod <- modele.fun(theta, obs) 
      fo <- forward(mod)
      ba <- backward(mod, fo)
      nb.fw <- nb.fw + 1
      nb.bw <- nb.bw + 1

      theta <- M.step.fun(obs, ba)
      U[,3] <- theta
      # *** keep trace of theta ***
      if(trace.theta) Theta[,k] <- theta
      k <- k+1;
      if(k > max.iter) {
        EXIT <- TRUE
        break
      }

      # check convergence
      if( sqrt(sum((U[,2] - U[,3])**2)) < epsilon ) {
        EXIT <- TRUE
        break
      }

      beta <- squarem.beta(U)
      if(beta <= 0) { # continue EM iterations
        U[,1:2] <- U[,2:3] # shift U
      } else break; # we got a beta > 0
    }
    if(EXIT) break

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
        nb.fw <- nb.fw + 1

        ll1 <- sum(log(colSums(fo$alpha * mod$p.emiss)))
      }
      if(ll1 >= ll0) { # accept proposal
        theta <- theta1

        # *** keep trace of theta ***
        if(trace.theta) Theta[,k] <- theta
        k <- k+1;
        if(k > max.iter) {
          EXIT <- TRUE
          break;
        }

        U[,1] <- theta # it's our new starting point

        # take advantage that the forward has been done already to finish the E step
        ba <- backward(mod, fo)
        nb.bw <- nb.bw + 1

        # M step
        theta <- M.step.fun(obs, ba)
        U[,2] <- theta
        # *** keep trace of theta ***
        if(trace.theta) Theta[,k] <- theta
        k <- k+1;
        if(k > max.iter) {
          EXIT <- TRUE
          break
        }

        # check convergence
        if( sqrt(sum((U[,1] - U[,2])**2)) < epsilon ) {
          EXIT <- TRUE
        }

        break # break from the backtracking loop, ready to compute U[,3]
      }
      # else, backtracking
      beta <- beta/2
      if(beta <= 0.25) { # we bactracked too far, going back to EM
        beta <- 0
        U[,1:2] <- U[,2:3] # shift U, readdy for going back to beginning of "big loop"
        break
      }
    }
    if(EXIT) break
  }
  # EXIT is true !
  R <- list(theta = theta, iter = k-1, forwards = nb.fw, backwards = nb.bw)
  if(trace.theta) R$Theta <- Theta[, 1:(k-1) ]
  R
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
