library(steveHMM)
library(salad)

path <- "/home/sidonie/Bureau/github/steveHMM/examples/"

source(paste0(path, "geyser_continu_stat.r"))

test.init.gc <- function(obs, it){
  qN <- matrix(nrow = it, ncol = 22)
  colnames(qN) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l", "tps")
  BW <- matrix(nrow = it, ncol = 22)
  colnames(BW) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l", "tps")
  SQ <- matrix(nrow = it, ncol = 24)
  colnames(SQ) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l","forward", "backward", "tps")
  QNEM <- matrix(nrow = it, ncol = 24)
  colnames(QNEM) <- c("iter", "a", "b", "muc", "mul", "muls", "sdc", "sdl", "sdls", "pi_c", "pi_l", "init_a", "init_b", "init_muc", "init_mul", "init_muls", "init_sdc", "init_sdl", "init_sdls", "init_pi_c", "init_pi_l", "forward", "backward", "tps")

  for(i in 1:it){
    print(i)
    set.seed(5*i+3)
    par <- c(runif(2), runif(3, 0, 6), runif(3,1,3))

    d <- Sys.time()
    tr <- try(qn <- quasi_newton(par, obs, modele.geyser.continu.stat, lower = rep(0.01,8), upper = c(0.99,0.99, rep(Inf,6))))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      qN[i,1] <- qn$counts[1]
      qN[i,2:9] <- qn$par
      qN[i,12:19] <- par
      qN[i,22] <- f-d
    }

    test.em <- function(par){
      em <- EM(par, obs, modele.geyser.continu.stat, M.step.geyser.continu.stat, max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(em)}
    d <- Sys.time()
    tr <- try(em <- test.em(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      BW[i,1] <- em$iter
      BW[i,2:9] <- em$theta
      BW[i,12:19] <- par
      BW[i,22] <- f-d
    }

    test.sq <- function(par){
      squarem <- SQUAREM(par, obs, modele.geyser.continu.stat, M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)), max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(squarem)}
    d <- Sys.time()
    tr <- try(squarem <- test.sq(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      SQ[i,1] <- squarem$iter
      SQ[i,2:9] <- squarem$theta
      SQ[i,12:19] <- par
      SQ[i,22] <- squarem$forwards
      SQ[i,23] <- squarem$backwards
      SQ[i,24] <- f-d
    }

    test.qnem <- function(par){
      qnem <- QNEM(par, obs, modele.geyser.continu.stat, M.step.geyser.continu.stat, lower = rep(0,8), upper = c(1,1, rep(Inf,6)) )
      return(qnem)}
    d <- Sys.time()
    tr <- try(qnem <- test.qnem(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      QNEM[i,1] <- qnem$iter
      QNEM[i,2:9] <- qnem$theta
      QNEM[i,12:19] <- par
      QNEM[i,22] <- qnem$forwards
      QNEM[i,23] <- qnem$backwards
      QNEM[i,24] <- f-d
    }
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN, QNEM = QNEM))
}


#observations
library(MASS)
X.geyser <- faithful$eruptions

#test du nombre d'iterations necessaires avec des points de départ différents
set.seed(28)
nb_it_gc <- test.init.gc(X.geyser, 1000)
# saveRDS(nb_it_gc, paste0(path,"nb_it_geyser_continu_qnem_meot.rds"))