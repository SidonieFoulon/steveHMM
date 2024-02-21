quasi_newton <- function(theta, obs, modele.fun, trace = FALSE, ...) {
  last_theta <- theta
  last_value <- neg_log_likelihood_gradient(theta, obs, modele.fun)

  ff <- function(theta) {
    if(trace) cat("theta =", theta, "\n")
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

capture_quasi_newton <- function(...) {

  L <- list(...)
  L <- c(L, list(control = list(trace = 6)))
  trace <- capture.output( do.call(quasi_newton, L) )

  tr1 <- trace[grep( "^(Cauchy X|theta) =", trace ) ]
  pts <- t(sapply(tr1, getval, USE.NAMES=FALSE))
  rownames(pts) <- sapply(tr1, getname)
  pts
}

getval <- function(x) { 
  x <- strsplit(x, split = '\\s+')[[1]] ;  
  x <- x[ -(1:which(x == "=")) ]; 
  as.numeric(x) 
}
getname <- function(x) { 
  x <- strsplit(x, split = '\\s+')[[1]] ;  
  x[1]
}

