#' SQUAREM algorithm.
#'
#'
#' @param theta the initial parameters
#' @param obs the observations data
#' @param modele.fun a model function
#' @param M.step.fun a M step function :
#'        - input : obs, backward
#'        - output : theta
#' @param lower the lower bound of the space of parameters
#' @param upper the upper bound of the space of parameters
#' @param max.iter maximum number of iteration of the algorithm (default is 100)
#' @param trace.theta whether you want to keep theta estimation for each iteration (default is TRUE)
#' @param epsilon if criteria = "eps", algorithm converges if the difference with the previous iteration is lower than eps (default is 1e-5)
#' @param reltol if criteria = "reltol", constant related to the stopping criterion, often depending on machine precision (default is sqrt(.Machine$double.eps))
#' @param criteria stopping criterion (default is "reltol")
#'
#' @return This function returns the final estimation of the parameters in "theta", the number of iterations in "iter", the number of forward and backward algorithms steps in resp. "forward" and "backward", and the final likelihood in "likelihood". If trace.theta = TRUE, it will also return the parameters estimated for each iteration in "Theta".
#'
#'
#' @export




SQUAREM <- function(theta, obs, modele.fun, M.step.fun, lower, upper, max.iter = 100, trace.theta = TRUE, epsilon = 1e-5, reltol = sqrt(.Machine$double.eps), criteria = c("reltol", "eps")) {
  if(is.infinite(max.iter)) trace.theta <- FALSE
  criteria <- match.arg(criteria)

  l <- length(obs)

  d <- length(theta)
  if(missing(lower)) lower <- rep(-Inf, d)
  if(missing(upper)) upper <- rep(Inf, d)
  if(length(upper) != d) stop("upper and theta should have same length")
  if(length(lower) != d) stop("lower and theta should have same length")
  if(any(theta < lower | theta > upper)) stop("Initial value not in the bounding box")

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
  if(any(is.infinite(mod$p.emiss))) {
    stop("Infinite density in model at starting value")
  }

  fo <- forward(mod)
  ll <- fo$likelihood
  ba <- backward(fo)
  nb.fw <- nb.fw + 1L
  nb.bw <- nb.bw + 1L

  theta1 <- M.step.fun(obs, ba)
  theta <- ifelse(is.na(theta1), theta, theta1)
  U[,2] <- theta
  # *** keep trace of theta ***
  if(trace.theta) Theta[,2] <- theta
  k <- 3

  EXIT <- FALSE
  repeat { # The big loop

    repeat { # iterate EM until beta > 0
      mod <- modele.fun(theta, obs)
      if(any(is.infinite(mod$p.emiss))) {
        warning("Infinite density in model")
        EXIT <- TRUE
        break
      }

      fo <- forward(mod)
      ll1 <- fo$likelihood
      rel.ll <- abs(ll - ll1) / (abs(ll) + reltol)
      ba <- backward(fo)
      if(any(is.na(ba$phi))) {
        warning("Backward step failed")
        EXIT <- TRUE
        break
      }

      nb.fw <- nb.fw + 1L
      nb.bw <- nb.bw + 1L

      theta1 <- M.step.fun(obs, ba)
      theta <- ifelse(is.na(theta1), theta, theta1)
      U[,3] <- theta
      # *** keep trace of theta ***
      if(trace.theta) Theta[,k] <- theta
      k <- k+1;
      ll <- ll1
      if(k > max.iter) {
        EXIT <- TRUE
        break
      }

      # check convergence
      if(criteria == "eps") {
        if( sqrt(sum((U[,2] - U[,3])**2)) < epsilon ) {
          EXIT <- TRUE
          break
        }
      } else {
        if(rel.ll < reltol){
          EXIT <- TRUE
          break
        }
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
        ll01 <- -Inf
      } else {
        mod <- modele.fun(theta1, obs)
        if(any(is.infinite(mod$p.emiss))) {
          warning("Infinite density in model")
          EXIT <- TRUE
          break
        }
        fo <- forward(mod)
        nb.fw <- nb.fw + 1L

        ll01 <- sum(log(colSums(fo$alpha * mod$p.emiss)))
        rel.ll0 <- abs(ll0 - ll01) / (abs(ll0) + reltol)
      }
      if(ll01 >= ll0) { # accept proposal
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
        ba <- backward(fo)
        if(any(is.na(ba$phi))) {
          warning("Backward step failed")
          EXIT <- TRUE
          break
        }

        nb.bw <- nb.bw + 1L

        # M step
        theta1 <- M.step.fun(obs, ba)
        theta <- ifelse(is.na(theta1), theta, theta1)
        U[,2] <- theta
        # *** keep trace of theta ***
        if(trace.theta) Theta[,k] <- theta
        k <- k+1;
        if(k > max.iter) {
          EXIT <- TRUE
          break
        }

        # check convergence
        if(criteria == "eps") {
          if( sqrt(sum((U[,1] - U[,2])**2)) < epsilon ) {
            EXIT <- TRUE
          }
        }
        if(criteria == "reltol" ) {
          if(rel.ll0 < reltol){
            EXIT <- TRUE
          }
        }

        break # break from the backtracking loop, ready to compute U[,3]
      }
      # else, backtracking
      beta <- beta/2
      if(beta <= 0.25) { # we bactracked too far, going back to EM
        beta <- 0
        U[,1:2] <- U[,2:3] # shift U, ready for going back to beginning of "big loop"
        break
      }
    }
    ll <- ll01
    if(EXIT) break
  }
  # EXIT is true !
  R <- list(theta = theta, iter = k-1, forwards = nb.fw, backwards = nb.bw, likelihood = ll)
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
