#' QNEM algorithm.
#'
#'
#' @param theta the initial parameters
#' @param obs the observations data
#' @param modele.fun a model function
#' @param M.step.fun a M step function :
#'        - input : obs, backward
#'        - output : theta
#' @param max.iter maximum number of iteration of the algorithm (default is 100)
#' @param lin.coeff linear coefficients for constraints
#' @param lin.upper upper constraints
#' @param trace.theta whether you want to keep theta estimation for each iteration (default is TRUE)
#' @param reltol if criteria = "reltol", constant related to the stopping criterion, often depending on machine precision (default is sqrt(.Machine$double.eps))
#' @param nb.em the number of Baum-Welch EM algorithm to perform before switching to BFGS algorithm
#' @param verbose if \code{TRUE}, displays information on the process (default is FALSE)
#'
#' @return This function returns the final estimation of the parameters in "theta", the final negative likelihood in "neg.ll", the number of iterations in "iter" and the number of forward and backward algorithms steps in resp. "forward" and "backward". If trace.theta = TRUE, it will also return the parameters estimated for each iteration in "Theta".
#'
#' @export


QNEM2 <- function(theta, obs, modele.fun, M.step.fun, max.iter = 100, lin.coeff, lin.upper, 
                 trace.theta = TRUE, reltol = sqrt(.Machine$double.eps), nb.em, verbose = FALSE){
  c.armijo <- 1e-4
  tau_backt <- .1
  if(missing(nb.em)) {
    if(verbose) cat("number of EM iterations determined by convexity\n")
    auto.em <- TRUE
    nb.em <- 0
  } else {
    if(verbose) cat("number of EM iterations determined by user\n")
    auto.em <- FALSE
  }

  d <- length(theta)

  if(!missing(lin.coeff)) {
    if(ncol(lin.coeff) != d) stop("lin.coeff and theta dimensions are not compatible")
    if(length(lin.upper) != nrow(lin.coeff)) stop("lin.coeff and lin.upper dimensions are not compatible")
  } else {
    lin.coeff <- matrix(0, nrow = 0, ncol = d)
    lin.upper <- numeric(0)
  }
  nc <- nrow(lin.coeff)

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

  # valeur initiale pour H
  # H <- diag(d)
  H <- matrix(NA_real_, d, d)

  # On amorce les calculs avant d'entrer dans la boucle
  # ** fo est supposé exister en début de boucle pour le point courant **
  # take modele derivatives
  mod <- modele_derivatives(modele.fun, theta, obs)
  fo <- forward_ll(mod, TRUE)
  nb.fw <- nb.fw + 1L
  ll <- fo$likelihood
  gradient <- fo$likelihood.gradient
  convex <- FALSE

  # the big loop
  repeat {
    if(nb.em > 0 | (auto.em & !convex)) { # EM
      if(verbose) cat("EM step\n")
      # Finish the E step
      ba <- backward(fo)
      nb.bw <- nb.bw + 1L

      if(any(is.na(ba$phi))) {
        warning("Forward / backward step failed")
        i <- i-1
        break
      }

      # Etape M
      theta1 <- M.step.fun(obs, ba)

      # enforce contraints 
      C.theta <- lin.coeff %*% theta
      for(r in seq_len(nc)) {
        if( C.theta[r] > lin.upper[r] ) { 
          # back to the border
          theta <- theta + ( lin.upper[r] - C.theta[r] ) * lin.coeff[r,] / sum( lin.coeff[r,]**2 )
        }
      }

      # in some cases, the M step produces NA's (e.g. conditionnal probabilities computed as 0/0)
      # this should be a correct way to deal with this
      # when a latent state has probability 0,
      # some of the parameters may be undefined: keep last value
      theta1 <- ifelse(is.na(theta1), theta, theta1)

      # on commence à préparer la boucle suivante
      mod <- modele_derivatives(modele.fun, theta1, obs)
      fo <- forward_ll(mod, TRUE)

      nb.fw <- nb.fw + 1L
      ll1 <- fo$likelihood
      gradient1 <- fo$likelihood.gradient

      # can happen near the border
      if(any(is.na(gradient1))) {
        gradient1[] <- 0
      }
      if(is.na(ll1)) {
        warning("Forward / backward step failed")
        # this will interrupt the loop
        ll1 <- ll
      }
      # record changes
      s <- theta1 - theta
      y <- gradient1 - gradient
      rel.ll <- abs(ll - ll1) / (abs(ll) + reltol)

      # update H if convex
      if(sum(s*y) > 0) {
        if(verbose) cat("locally convex region\n")
        convex <- TRUE
        H <- H_update(H, s, y)
      }

      # update current point
      theta <- theta1
      gradient  <- gradient1
      ll <- ll1

      # prêt...
      nb.em <- nb.em - 1
      QN <- FALSE
    } else { # QN
      if(verbose) cat("QN\n")
      if(any(is.na(H))) { # aucun pas d'EM n'a updaté H
        H <- diag(d)
      }

      # research direction
      p <- -(H %*% gradient)

      # Linear constraints.
      C.theta <- lin.coeff %*% theta
      for(r in seq_len(nc)) {
        # back to the border if necessary
        if( C.theta[r] > lin.upper[r] ) { 
          theta <- theta + ( lin.upper[r] - C.theta[r] ) * lin.coeff[r,] / sum( lin.coeff[r,]**2 )
        }
        # simply project research direction
        # are we on this border ?
        if( abs( C.theta[r] - lin.upper[r] ) < 1e-8 ) { 
           p <- p - sum(lin.coeff[r, ] * p) * lin.coeff[r,] / sum( lin.coeff[r,]**2 )     
        }
      }

      p.grad <- sum( (p * gradient) )
      if(p.grad > 0) { # pas une bonne direction... ne devrait pas arriver
        if(verbose) cat("Bad direction in BFGS\n")
        # on va zapper la line search et repartir sur l'EM
        lambda.i.max <- rep(0, nc) # ceci va produire lambda.max = 0
      } else {
        # preparing for backtracking
        # lambda.i.max <- ifelse(p > 0, (upper - theta)/p, ifelse(p < 0, (lower - theta)/p, Inf))
        C.p <- lin.coeff %*% p
        lambda.i.max <- ifelse(C.p > 0, (lin.upper - C.theta)/C.p, Inf)
      }

      lambda.max <- min(lambda.i.max, na.rm = TRUE) 
      if(lambda.max < 0) { # ne devrait pas arriver
        if(verbose) cat("BFGS failing...\n") 
        lambda.max <- 0
      }

      lambda <- min(1, lambda.max)
      repeat { # bracktracking loop
        theta1 <- theta + lambda * p
        if(all(theta1 == theta)) # did we really backtrack this far ? (or lambda was 0!)
          break
        if(lambda == lambda.max) { # block variables ???
        }
        mod <- modele_derivatives(modele.fun, theta1, obs)
        fo1 <- forward_ll(mod, TRUE)
        nb.fw <- nb.fw + 1L
        # likelihood undefined -> backtrack
        ll1 <- if(is.na(fo1$likelihood)) Inf else fo1$likelihood
        gradient1 <- fo1$likelihood.gradient
        rel.ll <- abs(ll - ll1) / (abs(ll) + reltol)
        # Armijo rule
        if(ll1 <= ll + c.armijo * lambda * p.grad)
          break
        lambda <- lambda * tau_backt
        if(verbose) cat("backtracking\n")
      }
      # out of the backtracking loop with a "good" theta1
      # update H
      s <- theta1 - theta
      y <- gradient1 - gradient
      if(sum(s*y) <= 0) { # curvature condition [catches also the case theta1 = theta]
        if(verbose) cat("Non convex spot\n")
        convex <- FALSE
        # nb.em <- 1
      } else {
        H <- H_update(H, s, y)
        # ------ if new wariables are blocked -----
      }
      # update loop variables
      ll <- ll1
      fo <- fo1
      theta <- theta1
      gradient <- fo$likelihood.gradient # "raw" gradient ! the beginning of the loop takes care of this
      QN <- TRUE
    }
    if(trace.theta) Theta[,i] <- theta
    if((QN & auto.em & !convex) | (QN & nb.em > 0)) { # in QN, decided to force at least one EM step
      H <- matrix(NA_real_, d, d)
    } else {
      if(rel.ll < reltol | i == max.iter) break
    }
    i <- i+1
    if(verbose) cat("theta =", as.vector(theta), "\n")
    if(verbose) cat("neg log likelihood =", ll, "\n")
  }
  R <- list(theta = theta, neg.ll = ll, iter = i, forwards = nb.fw, backwards = nb.bw)
  if(trace.theta) R$Theta <- Theta[, 1:i ]
  R
}

