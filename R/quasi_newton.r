quasi_newton <- function(theta, obs, modele.fun, trace = FALSE, ...) {
  last_theta <- theta
  last_value <- neg_log_likelihood_gradient(theta, obs, modele.fun)
 
  ff <- function(theta) {
    if(trace) { cat("theta =", theta, "\n") }
    if(!all(theta == last_theta)) {
      last_theta <<- theta
      last_value <<- neg_log_likelihood_gradient(theta, obs, modele.fun)
    }
    return(last_value$value)
  }

  dff <- function(theta) {
    if(!all(theta == last_theta)) {
      last_theta <<- theta
      last_value <<- neg_log_likelihood_gradient(theta, obs, modele.fun)
    }
    return(last_value$gradient)
  }
  optim(theta, ff, dff, method = "L-BFGS-B", ...)
}
