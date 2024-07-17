library(steveHMM)
library(salad)

#to adapt to your environnment
path <- "/home/sidonie/Bureau/github/steveHMM/examples/HBD_segments/"

source(paste0(path, "HBD_example.r"))

test.init.hbd <- function(obs, it){
  qN <- matrix(nrow = it, ncol = 7)
  colnames(qN) <- c("iter", "f", "a", "init_f", "init_a", "tps", "likelihood")
  BW <- matrix(nrow = it, ncol = 7)
  colnames(BW) <- c("iter", "f", "a", "init_f", "init_a", "tps", "likelihood")
  SQ <- matrix(nrow = it, ncol = 9)
  colnames(SQ) <- c("iter", "f", "a", "init_f", "init_a", "forward", "backward", "tps", "likelihood")
  QNEM <- matrix(nrow = it, ncol = 9)
  colnames(QNEM) <- c("iter", "f", "a", "init_f", "init_a", "forward", "backward", "tps", "likelihood")

  for(i in 1:it){
    print(i)
    RNGkind("Mer") ; set.seed(5*i+3)
    par <- c(runif(1,0,0.5), runif(1,0,10))

    d <- Sys.time()
    tr <- try(qn <- quasi_newton(par, obs, modele.HBD.stat, lower = rep(0.01,2), upper = c(0.99,Inf)))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      qN[i,1] <- qn$counts[1]
      qN[i,2:3] <- qn$par
      qN[i,4:5] <- par
      qN[i,6] <- f-d
      qN[i,7] <- qn$value
    }

    test.em <- function(par){
      em <- EM(par, obs, modele.HBD.stat, M.step.HBD.stat, max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(em)}
    d <- Sys.time()
    tr <- try(em <- test.em(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      BW[i,1] <- em$iter
      BW[i,2:3] <- em$theta
      BW[i,4:5] <- par
      BW[i,6] <- f-d
      BW[i,7] <- em$likelihood
    }

    test.sq <- function(par){
      squarem <- SQUAREM(par, obs, modele.HBD.stat, M.step.HBD.stat, lower = rep(0,2), upper = c(1,Inf), max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(squarem)}
    d <- Sys.time()
    tr <- try(squarem <- test.sq(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      SQ[i,1] <- squarem$iter
      SQ[i,2:3] <- squarem$theta
      SQ[i,4:5] <- squarem$theta
      SQ[i,6] <- squarem$forwards
      SQ[i,7] <- squarem$backwards
      SQ[i,8] <- f-d
      SQ[i,9] <- squarem$likelihood
    }

    test.qnem <- function(par){
      qnem <- QNEM(par, obs, modele.HBD.stat, M.step.HBD.stat, lower = rep(0,2), upper = c(1,Inf), max.iter = Inf, trace.theta = FALSE )
      return(qnem)}
    d <- Sys.time()
    tr <- try(qnem <- test.qnem(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      QNEM[i,1] <- qnem$iter
      QNEM[i,2:3] <- qnem$theta
      QNEM[i,4:5] <- par
      QNEM[i,6] <- qnem$forwards
      QNEM[i,7] <- qnem$backwards
      QNEM[i,8] <- f-d
      QNEM[i,9] <- qnem$neg.ll
    }
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN, QNEM = QNEM))
}


# observations
D <- read.table(paste0(path,"mozza.tsv"), header = TRUE)
RNGkind("Mer") ; set.seed(28) ; I <- as.vector(tapply(as.double(seq_along(D$dist)), D$dist, function(x) if(length(x) == 1) x else sample(x, 1)))
X.HBD <- D$X[I]
p.HBD <- D$p[I]

#comparison test with 1000 different starting points
RNGkind("Mer") ; set.seed(28)
nb_it_hbd <- test.init.hbd(X.HBD, 1000)
saveRDS(nb_it_hbd, "/home/sidonie/Bureau/results/HMM/res_hbd.rds") #path to adapt to your environment
