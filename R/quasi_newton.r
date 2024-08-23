#' Quasi-Newton L-BFGS-B algorithm.
#'
#'
#' @param theta the initial parameters
#' @param obs the observations data
#' @param modele.fun a model function
#' @param trace whether you want to keep theta estimation for each iteration (default is FALSE)
#' @param ... extra parameters for 'optim'
#'
#' @details You can also add optim function parameters.
#'
#' @return This function returns the same values than the optim function. Moreover, the number of forward algorithm steps can be found in "nb.fw".
#'
#' @export

quasi_newton <- function(theta, obs, modele.fun, trace = FALSE, ...) {
  last_theta <- theta
  last_value <- neg_log_likelihood_gradient(theta, obs, modele.fun)
  nb.fw <- 1

  ff <- function(theta) {
    if(trace) { cat("theta =", theta, "\n") }
    if(!all(theta == last_theta)) {
      last_theta <<- theta
      last_value <<- neg_log_likelihood_gradient(theta, obs, modele.fun)
      nb.fw <<- nb.fw + 1
    }
    return(last_value$likelihood)
  }

  dff <- function(theta) {
    if(!all(theta == last_theta)) {
      last_theta <<- theta
      last_value <<- neg_log_likelihood_gradient(theta, obs, modele.fun)
      nb.fw <<- nb.fw + 1
    }
    return(last_value$likelihood.gradient)
  }
  o <- optim(theta, ff, dff, method = "L-BFGS-B", ...)
  o$forwards <- nb.fw
  o
}
