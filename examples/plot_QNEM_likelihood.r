# Plot QNEM likelihood

#libraries
library(steveHMM)
library(salad)

###UMBRELLA###
#load data and functions
path <- "/home/sidonie/Bureau/github/steveHMM/examples/umbrella/"
source(paste0(path, "umbrella_example.r"))

#draw initial parameters
RNGkind("Mer") ; set.seed(23)
par <- runif(2)

#run QNEM
qnem.um <- QNEM(theta = par, obs = X.parapluie, modele.fun = modele.parapluie, M.step.fun = M.step.parapluie, lower = rep(0,2), upper = c(1,1), max.iter = Inf, trace.theta = TRUE, verbose = TRUE )

#get neg log likelihood from verbose
trace.um <- capture.output(QNEM(theta = par, obs = X.parapluie, modele.fun = modele.parapluie, M.step.fun = M.step.parapluie, lower = rep(0,2), upper = c(1,1), max.iter = Inf, trace.theta = TRUE, verbose = TRUE ))
trace.um1 <- trace.um[grep( "^(neg log likelihood) =", trace.um ) ]

getval <- function(x) {
  x <- strsplit(x, split = '\\s+')[[1]] ;
  x <- x[ -(1:which(x == "=")) ];
  as.numeric(x)
}

pts.um <- t(sapply(trace.um1, getval, USE.NAMES=FALSE))
pts.um <- as.vector(pts.um)


###GEYSER DICHO###
#load data and functions
path <- "/home/sidonie/Bureau/github/steveHMM/examples/geyser_dicho/"
source(paste0(path, "geyser_dicho_stat.r"))

#draw initial parameters
RNGkind("Mer") ; set.seed(23)
par <- runif(5)

#run QNEM
qnem.gd <- QNEM(theta = par, obs = X.dich, modele.fun = modele.geyser.stat, M.step.fun = M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = TRUE, verbose = TRUE)

#get neg log likelihood from verbose
trace.gd <- capture.output(QNEM(theta = par, obs = X.dich, modele.fun = modele.geyser.stat, M.step.fun = M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = TRUE, verbose = TRUE))
trace.gd1 <- trace.gd[grep( "^(neg log likelihood) =", trace.gd ) ]

pts.gd <- t(sapply(trace.gd1, getval, USE.NAMES=FALSE))
pts.gd <- as.vector(pts.gd)



###GEYSER CONT###
#load data and functions
path <- "/home/sidonie/Bureau/github/steveHMM/examples/geyser_continu/"
source(paste0(path, "geyser_continu_stat.r"))

#draw initial parameters
RNGkind("Mer") ; set.seed(23)
par <- c(runif(2), runif(3, 0, 6), runif(3,1,3))

#run QNEM
qnem.gc <- QNEM(theta = par, obs = X.geyser, modele.fun = modele.geyser.continu.stat, M.step.fun = M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = TRUE, verbose = TRUE )

#get neg log likelihood from verbose
trace.gc <- capture.output(QNEM(theta = par, obs = X.geyser, modele.fun = modele.geyser.continu.stat, M.step.fun = M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = TRUE, verbose = TRUE ))
trace.gc1 <- trace.gc[grep( "^(neg log likelihood) =", trace.gc ) ]

pts.gc <- t(sapply(trace.gc1, getval, USE.NAMES=FALSE))
pts.gc <- as.vector(pts.gc)



###HBD###
#load data and functions
path <- "/home/sidonie/Bureau/github/steveHMM/examples/HBD_segments/"
source(paste0(path, "HBD_example.r"))

#draw initial parameters
RNGkind("Mer") ; set.seed(23)
par <- c(runif(1,0,0.5), runif(1,0,10))

#run QNEM
qnem.hbd <- QNEM(theta = par, obs = X.HBD, modele.fun = modele.HBD.stat, M.step.fun = M.step.HBD.stat, lower = rep(0,2), upper = c(1,Inf), max.iter = Inf, trace.theta = TRUE, verbose = TRUE )

#get neg log likelihood from verbose
trace.hbd <- capture.output(QNEM(theta = par, obs = X.HBD, modele.fun = modele.HBD.stat, M.step.fun = M.step.HBD.stat, lower = rep(0,2), upper = c(1,Inf), max.iter = Inf, trace.theta = TRUE, verbose = TRUE ))
trace.hbd1 <- trace.hbd[grep( "^(neg log likelihood) =", trace.hbd ) ]

pts.hbd <- t(sapply(trace.hbd1, getval, USE.NAMES=FALSE))
pts.hbd <- as.vector(pts.hbd)


###ALL PLOTS###
par(mfrow = c(2,2))
plot(pts.um, type = "o", main = "(A) Umbrella example", xlab = "Iterations", ylab = "Negative log-likelihood")
plot(pts.gd, type = "o", main = "(B) Dichotomised geyser example", xlab = "Iterations", ylab = "Negative log-likelihood")
plot(pts.gc, type = "o", main = "(C) Continuous geyser example", xlab = "Iterations", ylab = "Negative log-likelihood")
plot(pts.hbd, type = "o", main = "(D) HBD segments example", xlab = "Iterations", ylab = "Negative log-likelihood")
