#' Baum-Welch EM algorithm.
#'
#'
#' @param theta the initial parameters
#' @param obs the observations data
#' @param modele.fun a model function
#' @param M.step.fun a M step function :
#'        - input : obs, backward
#'        - output : theta
#' @param max.iter maximum number of iteration of the algorithm (default is 100)
#' @param trace.theta whether you want to keep theta estimation for each iteration (default is TRUE)
#' @param epsilon if criteria = "eps", algorithm converges if the difference with the previous iteration is lower than eps (default is 1e-5)
#' @param reltol if criteria = "reltol", constant related to the stopping criterion, often depending on machine precision (default is sqrt(.Machine$double.eps))
#' @param criteria stopping criterion (default is "reltol")
#'
#' @return This function returns the final estimation of the parameters in "theta", the number of iterations in "iter" and the final likelihood in "likelihood". If trace.theta = TRUE, it will also return the parameters estimated for each iteration in "Theta".
#'
#'
#' @export


EM <- function(theta, obs, modele.fun, M.step.fun, max.iter = 100, trace.theta = TRUE, epsilon = 1e-5, reltol = sqrt(.Machine$double.eps), criteria = c("reltol", "eps")){
  if(is.infinite(max.iter)) trace.theta <- FALSE
  max.iter <- max(max.iter, 3)
  criteria <- match.arg(criteria)

  l <- length(obs)

  if(trace.theta) {
    Theta <- matrix(NA_real_, ncol = max.iter, nrow = length(theta))
    rownames(Theta) <- names(theta)
    Theta[,1] <- theta
  }

  i <- 2
  mod <- modele.fun(theta, obs)
  fo <- forward(mod)
  ll <- fo$likelihood
  ba <- backward(fo)
  theta <- M.step.fun(obs, ba)
  if(trace.theta) Theta[,i] <- theta

  i <- 3

  repeat {
    mod <- modele.fun(theta, obs)
    if(any(is.infinite(mod$p.emiss))) {
      warning("Infinite density in model")
      i <- i-1
      break
    }
    # Etape E
    fo <- forward(mod)
    ll1 <- fo$likelihood
    rel.ll <- abs(ll - ll1) / (abs(ll) + reltol)
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
    ll <- ll1

    if(trace.theta) Theta[,i] <- theta
    if(criteria == "eps") {
      if(e < epsilon | i == max.iter) break;
    } else {
      if(rel.ll < reltol | i == max.iter) break;
    }
    i <- i+1
  }
  R <- list(theta = theta, iter = i, likelihood = ll)
  if(trace.theta) R$Theta <- Theta[, 1:i ]
  R
}


