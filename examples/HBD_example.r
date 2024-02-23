## on note : S in {nHBD : 1 ; HBD : 2}
##           X in {AA : 1 ; Aa : 2 ; aa : 3}

modele.HBD <- function(theta, obs, name.S = c("nHBD", "HBD"), d = 0.1, pa = runif(length(obs)), eps = 0.01) { #theta takes f and a
  # transition matrix
  f <- theta[1]
  a <- theta[2]
  trans <- matrix( c((1-exp(-a*d))*(1-f)+exp(-a*d), (1-exp(-a*d))*f, (1-exp(-a*d))*(1-f), (1-exp(-a*d))*f+exp(-a*d)), byrow = TRUE, nrow = length(name.S), ncol = length(name.S))
  colnames(trans) <- rownames(trans) <- name.S

  # stationnary distribution
  x <- (f * (exp(-a*d) - 1) - exp(-a*d) + 1) / (1 - exp(-a*d))
  y <- 1-x
  stat <- c(nHBD = x, HBD = y)

  # emission probas are not re estimated : depends on alternative allele freq pa and genotypic error eps
  pA <- 1-pa
  p.emiss <- rbind(ifelse(obs == 1, pA**2, ifelse(obs == 2, 2*pA*pa, ifelse(obs == 3, pa**2, 1))),
                   ifelse(obs == 1, (1-eps)*pA + eps*(pA**2), ifelse(obs == 2, eps*2*pA*pa, ifelse(obs == 3, (1-eps)*pa + eps*(pa**2), 1))))
  rownames(p.emiss) <- name.S
  list(trans = trans, pi = stat, p.emiss = p.emiss)
}

M.step.HBD <- function(obs, backward, d=0.1) {
  l <- ncol(backward$phi)

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


#test on some random observations :

X.HBD <- sample(1:3, 100) #we expect f ~= 0 since it is random (no consiguinity)
theta.HBD <- c(0.05,0.05)

em.HBD <- EM(theta.HBD, X.HBD, modele.HBD, M.step.HBD, 10)
