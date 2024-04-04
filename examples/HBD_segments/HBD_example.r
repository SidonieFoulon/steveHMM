## on note : S in {nHBD : 1 ; HBD : 2}
##           X in {AA : 1 ; Aa : 2 ; aa : 3}

modele.HBD <- function(theta, obs, name.S = c("nHBD", "HBD"), d = 0.1, pa = p.HBD, eps = 0.01) { #theta takes f and a
  # transition matrix
  f <- theta[1]
  a <- theta[2]
  trans <- matrix( c((1-exp(-a*d))*(1-f)+exp(-a*d), (1-exp(-a*d))*f, (1-exp(-a*d))*(1-f), (1-exp(-a*d))*f+exp(-a*d)), byrow = TRUE, nrow = length(name.S), ncol = length(name.S))
  colnames(trans) <- rownames(trans) <- name.S

  # stationnary distribution
  nHBD <- theta[3]
  HBD <- theta[4]
  stat <- c(nHBD = 1-f, HBD = f)

  # emission probas are not re estimated : depends on alternative allele freq pa and genotypic error eps
  pA <- 1-pa
  p.emiss <- rbind(ifelse(obs == 1, pA**2, ifelse(obs == 2, 2*pA*pa, ifelse(obs == 3, pa**2, 1))),
                   ifelse(obs == 1, (1-eps)*pA + eps*(pA**2), ifelse(obs == 2, eps*2*pA*pa, ifelse(obs == 3, (1-eps)*pa + eps*(pa**2), 1))))
  rownames(p.emiss) <- name.S
  list(trans = trans, pi = stat, p.emiss = p.emiss)
}

M.step.HBD <- function(obs, backward, d=0.1) {
  l <- ncol(backward$phi)
  ba <- backward
  # first re estimate transition matrix
  t11 <- sum(ba$delta[1,1,-1]) / sum(ba$phi[1,-l])
  t12 <- sum(ba$delta[1,2,-1]) / sum(ba$phi[1,-l])
  t21 <- sum(ba$delta[2,1,-1]) / sum(ba$phi[2,-l])
  t22 <- sum(ba$delta[2,2,-1]) / sum(ba$phi[2,-l])
  tr <- matrix(c(t11, t12, t21, t22),byrow=TRUE, nrow=2)
  # f (-> trans)
  f <- tr[1,2] / (tr[1,2] + tr[2,1])
  # a (-> trans)
  a <- log(1 + tr[1,2] - tr[2,1]) / -d

  #no emission parameter to estimate

  #initial states to estimate
  nHBD <- backward$phi[1,1]
  HBD <- backward$phi[2,1]

  c(f, a, nHBD, HBD)
}

modele.HBD.stat <- function(theta, obs, name.S = c("nHBD", "HBD"), d = 0.1, pa = p.HBD, eps = 0.01) { #theta takes f and a
  # transition matrix
  f <- theta[1]
  a <- theta[2]
  trans <- matrix( c((1-exp(-a*d))*(1-f)+exp(-a*d), (1-exp(-a*d))*f, (1-exp(-a*d))*(1-f), (1-exp(-a*d))*f+exp(-a*d)), byrow = TRUE, nrow = length(name.S), ncol = length(name.S))
  colnames(trans) <- rownames(trans) <- name.S

  # stationnary distribution
  stat <- c(nHBD = 1-f, HBD = f)

  # emission probas are not re estimated : depends on alternative allele freq pa and genotypic error eps
  pA <- 1-pa
  p.emiss <- rbind(ifelse(obs == 1, pA**2, ifelse(obs == 2, 2*pA*pa, ifelse(obs == 3, pa**2, 1))),
                   ifelse(obs == 1, (1-eps)*pA + eps*(pA**2), ifelse(obs == 2, eps*2*pA*pa, ifelse(obs == 3, (1-eps)*pa + eps*(pa**2), 1))))
  rownames(p.emiss) <- name.S
  list(trans = trans, pi = stat, p.emiss = p.emiss)
}

M.step.HBD.stat <- function(obs, backward, d=0.1) {
  l <- ncol(backward$phi)
  ba <- backward
  # first re estimate transition matrix
  t11 <- sum(ba$delta[1,1,-1]) / sum(ba$phi[1,-l])
  t12 <- sum(ba$delta[1,2,-1]) / sum(ba$phi[1,-l])
  t21 <- sum(ba$delta[2,1,-1]) / sum(ba$phi[2,-l])
  t22 <- sum(ba$delta[2,2,-1]) / sum(ba$phi[2,-l])
  tr <- matrix(c(t11, t12, t21, t22),byrow=TRUE, nrow=2)
  # f (-> trans)
  f <- tr[1,2] / (tr[1,2] + tr[2,1])
  # a (-> trans)
  a <- log(1 + tr[1,2] - tr[2,1]) / -d

  #no emission parameter to estimate

  c(f,a)
}


#our observations :
library(Mozza)
if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
H <- read.bed.matrix(filepath)
set.seed(28)
x <- make.inbreds(1, H, f = 0.0625, a = 0.064, segments = TRUE)
x$inbred.coef # f = 0.078 on the whole genome, not the final value to estimate since we will not work on whole genome
x$segments # 2 HBD segments on chr 15
bm <- x$bed.matrix
# on récupère les fréquences dans la pop européenne
bm@p <- H@p
# SNPs fréquents du chr 17
bm <- select.snps(bm, chr == 17 & bm@p > 0.05 & bm@p < 0.95)
# un SNP tous les 0.1 cM environ
bm <- bm[, which(diff(round(bm@snps$dist,1))>0)]

X.HBD <- as.vector(as.matrix(bm))+1
p.HBD <- bm@p

par.HBD <- c(f = 0.05, a = 0.05, nHBD = 0.95, HBD = 0.05)

# qN
qn.hbd <- quasi_newton(par.HBD, X.HBD, modele.HBD, lower = c(.01,.01), upper = c(.99,Inf))
qn.hbd.trace <- capture_quasi_newton(par.HBD, X.HBD, modele.HBD, lower = c(.01,.01), upper = c(.99,Inf))

# EM
em.hbd <- EM(par.HBD, X.HBD, modele.HBD, M.step.HBD, trace.theta = TRUE)

# SQUAREM
squarem.hbd <- SQUAREM(par.HBD, X.HBD, modele.HBD, M.step.HBD, trace.theta = TRUE)

# QNEM
qnem.hbd <- QNEM(par.HBD, X.HBD, modele.HBD, M.step.HBD,trace.theta = TRUE, verbose = TRUE)


