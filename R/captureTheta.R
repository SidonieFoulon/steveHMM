captureTheta <- function(theta, optim.fun, method = "L-BFGS-B", low=rep(0,length(theta)) + 0.01, up=rep(1,length(theta)) - 0.01){

  trace <- capture.output(optim(par = theta, fn = optim.fun, method = method, lower=low, upper=up, control=list(trace=6)))

  trace <- trace[grep("X = ", trace)]

  trace <- trace[-grep("Cauchy", trace)]

  trace <- unlist(strsplit(trace, split = '\\s'))

  trace <- trace[-which(trace=="X")]

  trace <- trace[-which(trace=="=")]

  trace <- matrix(trace, nrow = length(theta), byrow = FALSE)

  trace <- cbind(theta, trace) # ajout de l'initialisation

  rownames(trace) <- names(theta)

  class(trace) <- "numeric"

  return(trace)

}


captureThetaCauchy <- function(theta, optim.fun, method = "L-BFGS-B", low=rep(0,length(theta)) + 0.01, up=rep(1,length(theta)) - 0.01){

  trace <- capture.output(optim(par = theta, fn = optim.fun, method = method, lower=low, upper=up, control=list(trace=6)))

  trace <- trace[grep("X = ", trace)]

  trace <- unlist(strsplit(trace, split = '\\s'))

  trace <- trace[-which(trace=="Cauchy")]

  trace <- trace[-which(trace=="X")]

  trace <- trace[-which(trace=="=")]

  trace <- trace[-which(trace=="")]

  trace <- matrix(trace, nrow = length(theta), byrow = FALSE)

  trace <- cbind(theta, trace) # ajout de l'initialisation

  rownames(trace) <- names(theta)

  class(trace) <- "numeric"

  return(trace)

}
