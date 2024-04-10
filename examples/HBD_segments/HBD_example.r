## on note : S in {nHBD : 1 ; HBD : 2}
##           X in {AA : 1 ; Aa : 2 ; aa : 3}

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
  list(trans = trans, pi = stat, p.emiss = p.emiss, theta = theta)
}

M.step.HBD.stat <- function(obs, backward, d=0.1) {
  l <- ncol(backward$phi)

  # first re estimate transition matrix
  a0 <- backward$trans[1,2]
  b0 <- backward$trans[2,1]
  D11 <- sum(backward$delta[1,1,])
  D12 <- sum(backward$delta[1,2,]) + backward$phi[2,1]
  D21 <- sum(backward$delta[2,1,]) + backward$phi[1,1]
  D22 <- sum(backward$delta[2,2,])

  repeat {
    a <- (D12 - a0*(1-a0)/(a0+b0)) / (D11 + D12)
    b <- (D21 - b0*(1-b0)/(a0+b0)) / (D21 + D22)
    if(a >= 0 & a <= 1 & b >= 0 & b <= 1) break
  }

  tr <- matrix(c(1 - a, a, b, 1 - b),byrow=TRUE, nrow=2)

  # f 
  f <- tr[1,2] / (tr[1,2] + tr[2,1])
  # a
  a <- log(1 - tr[1,2] - tr[2,1]) / -d

  #no emission parameter to estimate

  c(f,a)
}


#our observations :
D <- read.table("mozza.tsv", header = TRUE)

# a submap :
set.seed(28)
I <- as.vector(tapply(as.double(seq_along(D$dist)), D$dist, function(x) if(length(x) == 1) x else sample(x, 1)))
X.HBD <- D$X[I]
p.HBD <- D$p[I]

par.HBD <- c(f = 0.05, a = 0.05)

em <- EM(par.HBD, X.HBD, modele.HBD.stat, M.step.HBD.stat)
qnem <- QNEM(par.HBD, X.HBD, modele.HBD.stat, M.step.HBD.stat, upper = c(1,Inf), lower = c(0,0))
squarem <- SQUAREM(par.HBD, X.HBD, modele.HBD.stat, M.step.HBD.stat, upper = c(1,Inf), lower = c(0,0))
qn <- quasi_newton(par.HBD, X.HBD, modele.HBD.stat, upper = c(1,Inf), lower = c(0,0))

