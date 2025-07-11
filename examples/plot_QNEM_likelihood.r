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


#run QN
qn.um <- quasi_newton(par, X.parapluie, modele.parapluie, lower = rep(0.01,2), upper = c(0.99,0.99), control=list(trace = 1, REPORT=1))

#get neg log likelihood from verbose
tr.qn.um <- capture.output(quasi_newton(par, X.parapluie, modele.parapluie, lower = rep(0.01,2), upper = c(0.99,0.99), control=list(trace = 1, REPORT=1)))
tr.qn.um1 <- tr.qn.um[grep( "^iter", tr.qn.um ) ]
tr.qn.um2<- strsplit(tr.qn.um1, split = '\\s+')
nll.qn.um <- vector()
for(i in 1:length(tr.qn.um2)) {nll.qn.um <- c(nll.qn.um, tr.qn.um2[[i]][4])}
nll.qn.um <- as.numeric(nll.qn.um)



#run SQUAREM
sq.um <- SQUAREM(par, X.parapluie, modele.parapluie, M.step.parapluie, lower = rep(0,2), upper = c(1,1), max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.sq.um <- vector()
for (i in 1:ncol(sq.um$Theta)){
  colonne <- sq.um$Theta[,i]
  mod <- modele.parapluie(colonne, X.parapluie)
  fo <- forward(mod)
  nll.sq.um <- c(nll.sq.um, -fo$likelihood)
}



#run EM
em.um <- EM(par, X.parapluie, modele.parapluie, M.step.parapluie, max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.em.um <- vector()
for (i in 1:ncol(em.um$Theta)){
  colonne <- em.um$Theta[,i]
  mod <- modele.parapluie(colonne, X.parapluie)
  fo <- forward(mod)
  nll.em.um <- c(nll.em.um, -fo$likelihood)
}





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



#run QN
qn.hbd <- quasi_newton(par, X.HBD, modele.HBD.stat, lower = rep(0.01,2), upper = c(0.99,Inf), control=list(trace = 1, REPORT=1))

#get neg log likelihood from verbose
tr.qn.hbd <- capture.output(quasi_newton(par, X.HBD, modele.HBD.stat, lower = rep(0.01,2), upper = c(0.99,Inf), control=list(trace = 1, REPORT=1)))
tr.qn.hbd1 <- tr.qn.hbd[grep( "^iter", tr.qn.hbd ) ]
tr.qn.hbd2<- strsplit(tr.qn.hbd1, split = '\\s+')
nll.qn.hbd <- vector()
for(i in 1:length(tr.qn.hbd2)) {nll.qn.hbd <- c(nll.qn.hbd, tr.qn.hbd2[[i]][4])}
nll.qn.hbd <- as.numeric(nll.qn.hbd)


#run SQUAREM
sq.HBD <- SQUAREM(par, X.HBD, modele.HBD.stat, M.step.HBD.stat, lower = rep(0,2), upper = c(1,Inf), max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.sq.hbd <- vector()
for (i in 1:ncol(sq.HBD$Theta)){
  colonne <- sq.HBD$Theta[,i]
  mod <- modele.HBD.stat(colonne, X.HBD)
  fo <- forward(mod)
  nll.sq.hbd <- c(nll.sq.hbd, -fo$likelihood)
}

#run EM
em.HBD <- EM(par, X.HBD, modele.HBD.stat, M.step.HBD.stat, max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.em.hbd <- vector()
for (i in 1:ncol(em.HBD$Theta)){
  colonne <- em$Theta[,i]
  mod <- modele.HBD.stat(colonne, X.HBD)
  fo <- forward(mod)
  nll.em.hbd <- c(nll.em.hbd, -fo$likelihood)
}



###GEYSER CONT###
#load data and functions
path <- "/home/sidonie/Bureau/github/steveHMM/examples/geyser_continu/"
source(paste0(path, "geyser_continu_stat.r"))

#draw initial parameters
RNGkind("Mer") ; set.seed(38)
par <- c(runif(2), runif(3, 0, 6), runif(3,1,3))

#run QNEM
qnem.gc <- QNEM(theta = par, obs = X.geyser, modele.fun = modele.geyser.continu.stat, M.step.fun = M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = TRUE, verbose = TRUE )

#get neg log likelihood from verbose
trace.gc <- capture.output(QNEM(theta = par, obs = X.geyser, modele.fun = modele.geyser.continu.stat, M.step.fun = M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = TRUE, verbose = TRUE ))
trace.gc1 <- trace.gc[grep( "^(neg log likelihood) =", trace.gc ) ]

pts.gc <- t(sapply(trace.gc1, getval, USE.NAMES=FALSE))
pts.gc <- as.vector(pts.gc)

#run QN
qn.gc <- quasi_newton(par, X.geyser, modele.geyser.continu.stat, lower = rep(0.01,8), upper = c(0.99,0.99, rep(Inf,6)), control=list(trace = 1, REPORT=1))

#get neg log likelihood from verbose
tr.qn.gc <- capture.output(quasi_newton(par, X.geyser, modele.geyser.continu.stat, lower = rep(0.01,8), upper = c(0.99,0.99, rep(Inf,6)), control=list(trace = 1, REPORT=1)))
tr.qn.gc1 <- tr.qn.gc[grep( "^iter", tr.qn.gc ) ]
tr.qn.gc2<- strsplit(tr.qn.gc1, split = '\\s+')
nll.qn.gc <- vector()
for(i in 1:length(tr.qn.gc2)) {nll.qn.gc <- c(nll.qn.gc, tr.qn.gc2[[i]][4])}
nll.qn.gc <- as.numeric(nll.qn.gc)

#run SQUAREM
sq.gc <- SQUAREM(par, X.geyser, modele.geyser.continu.stat, M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.sq.gc <- vector()
for (i in 1:ncol(sq.gc$Theta)){
  colonne <- sq.gc$Theta[,i]
  mod <- modele.geyser.continu.stat(colonne, X.geyser)
  fo <- forward(mod)
  nll.sq.gc <- c(nll.sq.gc, -fo$likelihood)
}

#run EM
em.gc <- EM(par, X.geyser, modele.geyser.continu.stat, M.step.geyser.continu.stat, max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.em.gc <- vector()
for (i in 1:ncol(em.gc$Theta)){
  colonne <- em.gc$Theta[,i]
  mod <- modele.geyser.continu.stat(colonne, X.geyser)
  fo <- forward(mod)
  nll.em.gc <- c(nll.em.gc, -fo$likelihood)
}




###GEYSER DICHO###
#load data and functions
path <- "/home/sidonie/Bureau/github/steveHMM/examples/geyser_dicho/"
source(paste0(path, "geyser_dicho_stat.r"))

#draw initial parameters
RNGkind("Mer") ; set.seed(43)
par <- runif(5)

#run QNEM
qnem.gd <- QNEM(theta = par, obs = X.dich, modele.fun = modele.geyser.stat, M.step.fun = M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = TRUE, verbose = TRUE)

#get neg log likelihood from verbose
trace.gd <- capture.output(QNEM(theta = par, obs = X.dich, modele.fun = modele.geyser.stat, M.step.fun = M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = TRUE, verbose = TRUE))
trace.gd1 <- trace.gd[grep( "^(neg log likelihood) =", trace.gd ) ]

pts.gd <- t(sapply(trace.gd1, getval, USE.NAMES=FALSE))
pts.gd <- as.vector(pts.gd)

#run QN
qn.gd <- quasi_newton(par, X.dich, modele.geyser.stat, lower = rep(0.01,5), upper = rep(0.99,5), control=list(trace = 1, REPORT=1))

#get neg log likelihood from verbose
tr.qn.gd <- capture.output(quasi_newton(par, X.dich, modele.geyser.stat, lower = rep(0.01,5), upper = rep(0.99,5), control=list(trace = 1, REPORT=1)))
tr.qn.gd1 <- tr.qn.gd[grep( "^iter", tr.qn.gd ) ]
tr.qn.gd2<- strsplit(tr.qn.gd1, split = '\\s+')
nll.qn.gd <- vector()
for(i in 1:length(tr.qn.gd2)) {nll.qn.gd <- c(nll.qn.gd, tr.qn.gd2[[i]][4])}
nll.qn.gd <- as.numeric(nll.qn.gd)

#run SQUAREM
sq.gd <- SQUAREM(par, X.dich, modele.geyser.stat, M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.sq.gd <- vector()
for (i in 1:ncol(sq.gd$Theta)){
  colonne <- sq.gd$Theta[,i]
  mod <- modele.geyser.stat(colonne, X.dich)
  fo <- forward(mod)
  nll.sq.gd <- c(nll.sq.gd, -fo$likelihood)
}

#run EM
em.gd <- EM(par, X.dich, modele.geyser.stat, M.step.geyser.stat, max.iter = 100, trace.theta = TRUE, criteria = "rel")
nll.em.gd <- vector()
for (i in 1:ncol(em.gd$Theta)){
  colonne <- em.gd$Theta[,i]
  mod <- modele.geyser.stat(colonne, X.dich)
  fo <- forward(mod)
  nll.em.gd <- c(nll.em.gd, -fo$likelihood)
}

###ALL PLOTS###
pdf(file = "/home/sidonie/Bureau/github/steveHMM/examples/likelihoods.Rplot.pdf", width = 11.69, height = 8.27)
par(mfrow = c(2,2))
plot(pts.um, type = "o", main = "(A) Umbrella example", xlab = "Iterations", ylab = "Negative log-likelihood", xlim = c(0,15))
lines(nll.sq.um, type = "o", col = "green")
lines(nll.em.um, type = "o", col = "purple")
lines(nll.qn.um, type = "o", col = "red")
plot(pts.gd, type = "o", main = "(B) Dichotomised geyser example", xlab = "Iterations", ylab = "Negative log-likelihood", xlim=c(0,25))
lines(nll.sq.gd, type = "o", col = "green")
lines(nll.em.gd, type = "o", col = "purple")
lines(nll.qn.gd, type = "o", col = "red")
legend("topright", pch = 1, lty = 1, legend = c("qN", "BW", "SQUAREM", "QNEM"), col = c("red", "purple", "green", "black"))
plot(pts.gc, type = "o", main = "(C) Continuous geyser example", xlab = "Iterations", ylab = "Negative log-likelihood", xlim = c(0,35))
lines(nll.sq.gc, type = "o", col = "green")
lines(nll.em.gc, type = "o", col = "purple")
lines(nll.qn.gc, type = "o", col = "red")
plot(pts.hbd, type = "o", main = "(D) HBD segments example", xlab = "Iterations", ylab = "Negative log-likelihood", xlim = c(0,45))
lines(nll.sq.hbd, type = "o", col = "green")
lines(nll.em.hbd, type = "o", col = "purple")
lines(nll.qn.hbd, type = "o", col = "red")
dev.off()
