library(steveHMM)
library(salad)

#to adapt to your environment
path <- "/home/sidonie/Bureau/github/steveHMM/"

source(paste0(path, "examples/umbrella/umbrella_example.r"))

test.init.u <- function(obs, it){
  qN <- matrix(nrow = it, ncol = 7)
  colnames(qN) <- c("iter", "a", "b", "init_a", "init_b", "tps", "likelihood")
  BW <- matrix(nrow = it, ncol = 7)
  colnames(BW) <- c("iter", "a", "b", "init_a", "init_b", "tps", "likelihood")
  SQ <- matrix(nrow = it, ncol = 9)
  colnames(SQ) <- c("iter", "a", "b", "init_a", "init_b", "forward", "backward", "tps", "likelihood")
  QNEM <- matrix(nrow = it, ncol = 9)
  colnames(QNEM) <- c("iter", "a", "b", "init_a", "init_b", "forward", "backward", "tps", "likelihood")

  for(i in 1:it){
    print(i)
    RNGkind("Mer") ; set.seed(5*i+3)
    par <- runif(2)

    d <- Sys.time()
    tr <- try(qn <- quasi_newton(par, obs, modele.parapluie, lower = rep(0.01,2), upper = c(0.99,0.99)))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      qN[i,1] <- qn$counts[1]
      qN[i,2:3] <- qn$par
      qN[i,4:5] <- par
      qN[i,6] <- f-d
      qN[i,7] <- qn$value
    }

    test.em <- function(par){
      em <- EM(par, obs, modele.parapluie, M.step.parapluie, max.iter = Inf, trace.theta = FALSE, criteria = "rel")
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
      squarem <- SQUAREM(par, obs, modele.parapluie, M.step.parapluie, lower = rep(0,2), upper = c(1,1), max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(squarem)}
    d <- Sys.time()
    tr <- try(squarem <- test.sq(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      SQ[i,1] <- squarem$iter
      SQ[i,2:3] <- squarem$theta
      SQ[i,4:5] <- par
      SQ[i,6] <- squarem$forwards
      SQ[i,7] <- squarem$backwards
      SQ[i,8] <- f-d
      SQ[i,9] <- squarem$likelihood
    }

    test.qnem <- function(par){
      qnem <- QNEM(par, obs, modele.parapluie, M.step.parapluie, lower = rep(0,2), upper = c(1,1), max.iter = Inf, trace.theta = FALSE )
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


#observations
X.parapluie <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1,
                 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2)

#test du nombre d'iterations necessaires avec des points de départ différents
RNGkind("Mer") ; set.seed(28)
nb_it_u <- test.init.u(X.parapluie, 1000)
#saveRDS(nb_it_u, paste0(path,"results/nb_it_umbrella_qnem_meot.rds"))
