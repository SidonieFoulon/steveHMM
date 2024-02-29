capture_quasi_newton <- function(...) {

  L <- list(...)
  L <- c(L, list(control = list(trace = 6)))
  trace <- capture.output( do.call(quasi_newton, L) )

  tr1 <- trace[grep( "^(Cauchy X|theta) =", trace ) ]
  pts <- sapply(tr1, getval, USE.NAMES=FALSE)
  colnames(pts) <- sapply(tr1, getname)
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

capture_optim <- function(...) {

  L <- list(...)
  L <- c(L, list(control = list(trace = 6)))
  trace <- capture.output( do.call(optim, L) )

  tr1 <- trace[grep( "^(Cauchy X|theta) =", trace ) ]
  pts <- t(sapply(tr1, getval, USE.NAMES=FALSE))
  rownames(pts) <- sapply(tr1, getname)
  pts
}
