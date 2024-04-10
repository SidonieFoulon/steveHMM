library(steveHMM)
library(salad)

#to adapt to your environment
path <- "/home/sidonie/Bureau/github/steveHMM/"

source(paste0(path, "examples/geyser_dicho/geyser_dicho_stat.r"))

test.init.gd <- function(obs, it){
  qN <- matrix(nrow = it, ncol = 17)
  colnames(qN) <- c("iter", "a", "b", "c", "d", "e", "pi_c", "pi_l", "init_a", "init_b", "init_c", "init_d", "init_e", "init_pi_c", "init_pi_l", "tps", "likelihood")
  BW <- matrix(nrow = it, ncol = 17)
  colnames(BW) <- c("iter", "a", "b", "c", "d", "e", "pi_c", "pi_l", "init_a", "init_b", "init_c", "init_d", "init_e", "init_pi_c", "init_pi_l", "tps", "likelihood")
  SQ <- matrix(nrow = it, ncol = 19)
  colnames(SQ) <- c("iter", "a", "b", "c", "d", "e", "pi_c", "pi_l", "init_a", "init_b", "init_c", "init_d", "init_e", "init_pi_c", "init_pi_l", "forward", "backward", "tps", "likelihood")
  QNEM <- matrix(nrow = it, ncol = 19)
  colnames(QNEM) <- c("iter", "a", "b", "c", "d", "e", "pi_c", "pi_l", "init_a", "init_b", "init_c", "init_d", "init_e", "init_pi_c", "init_pi_l", "forward", "backward", "tps", "likelihood")

  for(i in 1:it){
    print(i)
    RNGkind("Mer") ; set.seed(5*i+3)
    par <- runif(5)

    d <- Sys.time()
    tr <- try(qn <- quasi_newton(par, obs, modele.geyser.stat, lower = rep(0.01,5), upper = rep(0.99,5)))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      qN[i,1] <- qn$counts[1]
      qN[i,2:6] <- qn$par
      qN[i,9:13] <- par
      qN[i,16] <- f-d
      qN[i,17] <- qn$value
    }

    test.em <- function(par){
      em <- EM(par, obs, modele.geyser.stat, M.step.geyser.stat, max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(em)}
    d <- Sys.time()
    tr <- try(em <- test.em(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      BW[i,1] <- em$iter
      BW[i,2:6] <- em$theta
      BW[i,9:13] <- par
      BW[i,16] <- f-d
      BW[i,17] <- em$likelihood
    }

    test.sq <- function(par){
      squarem <- SQUAREM(par, obs, modele.geyser.stat, M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(squarem)}
    d <- Sys.time()
    tr <- try(squarem <- test.sq(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      SQ[i,1] <- squarem$iter
      SQ[i,2:6] <- squarem$theta
      SQ[i,9:13] <- par
      SQ[i,16] <- squarem$forwards
      SQ[i,17] <- squarem$backwards
      SQ[i,18] <- f-d
      SQ[i,19] <- squarem$likelihood
    }

    test.qnem <- function(par){
      qnem <- QNEM(par, obs, modele.geyser.stat, M.step.geyser.stat, lower = rep(0,5), upper = rep(1,5), max.iter = Inf, trace.theta = FALSE)
      return(qnem)}
    d <- Sys.time()
    tr <- try(qnem <- test.qnem(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      QNEM[i,1] <- qnem$iter
      QNEM[i,2:6] <- qnem$theta
      QNEM[i,9:13] <- par
      QNEM[i,16] <- qnem$forwards
      QNEM[i,17] <- qnem$backwards
      QNEM[i,18] <- f-d
      QNEM[i,19] <- qnem$neg.ll
    }
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN, QNEM = QNEM))
}


#observations
library(MASS)
X.dicho <- ifelse(faithful$eruptions < 3, 1,2)

#test du nombre d'iterations necessaires avec des points de départ différents
RNGkind("Mer") ; set.seed(28)
nb_it_gd <- test.init.gd(X.dicho, 1000)
#saveRDS(nb_it_gd, paste0(path,"results/nb_it_geyser_dicho_qnem_meot.rds"))
