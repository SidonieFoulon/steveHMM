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
