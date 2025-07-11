library(steveHMM)
library(salad)

#to adapt to your environnment
path <- "/home/sidonie/Bureau/github/steveHMM/examples/eye_movements/"

source(paste0(path, "eye_movement_trans.R"))

test.eye <- function(X, par, pairs){
  it <- length(pairs)
  qN <- matrix(nrow = it, ncol = 31)
  colnames(qN) <- c("iter",
                    "t12", "t13", "t21", "t23", "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2",
                    "init_t12", "init_t13", "init_t21", "init_t23", "init_t31", "init_t32",
                    "init_e12", "init_e13", "init_e21", "init_e23", "init_e31", "init_e32",
                    "init_pi1", "init_pi2",
                    "tps", "likelihood")
  BW <- matrix(nrow = it, ncol = 31)
  colnames(BW) <- c("iter",
                    "t12", "t13", "t21", "t23", "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2",
                    "init_t12", "init_t13", "init_t21", "init_t23", "init_t31", "init_t32",
                    "init_e12", "init_e13", "init_e21", "init_e23", "init_e31", "init_e32",
                    "init_pi1", "init_pi2",
                    "tps", "likelihood")
  SQ <- matrix(nrow = it, ncol = 33)
  colnames(SQ) <- c("iter",
                    "t12", "t13", "t21", "t23", "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2",
                    "init_t12", "init_t13", "init_t21", "init_t23", "init_t31", "init_t32",
                    "init_e12", "init_e13", "init_e21", "init_e23", "init_e31", "init_e32",
                    "init_pi1", "init_pi2",
                    "forward", "backward", "tps", "likelihood")
  QNEM <- matrix(nrow = it, ncol = 33)
  colnames(QNEM) <- c("iter",
                      "t12", "t13", "t21", "t23", "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2",
                      "init_t12", "init_t13", "init_t21", "init_t23", "init_t31", "init_t32",
                      "init_e12", "init_e13", "init_e21", "init_e23", "init_e31", "init_e32",
                      "init_pi1", "init_pi2",
                      "forward", "backward", "tps", "likelihood")

  for(i in 1:it){
    print(pairs[i])
    #RNGkind("Mer") ; set.seed(5*i+3)
    obs <- subset(X, X$PAIR == pairs[i])$READMODE


    d <- Sys.time()
    tr <- try(qn <- quasi_newton(par, obs, modele.eye, lower = rep(0.01,14), upper = rep(0.99,14)))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      qN[i,1] <- qn$counts[1]
      qN[i,2:15] <- qn$par
      qN[i,16:29] <- par
      qN[i,30] <- f-d
      qN[i,31] <- qn$value
    }

    test.em <- function(par){
      em <- EM(par, obs, modele.eye, M.step.eye, max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(em)}
    d <- Sys.time()
    tr <- try(em <- test.em(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      BW[i,1] <- em$iter
      BW[i,2:15] <- em$theta
      BW[i,16:29] <- par
      BW[i,30] <- f-d
      BW[i,31] <- em$likelihood
    }

    test.sq <- function(par){
      squarem <- SQUAREM(par, obs, modele.eye, M.step.eye, lower = rep(0,14), upper = rep(1,14), max.iter = Inf, trace.theta = FALSE, criteria = "rel")
      return(squarem)}
    d <- Sys.time()
    tr <- try(squarem <- test.sq(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      SQ[i,1] <- squarem$iter
      SQ[i,2:15] <- squarem$theta
      SQ[i,16:29] <- par
      SQ[i,30] <- squarem$forwards
      SQ[i,31] <- squarem$backwards
      SQ[i,32] <- f-d
      SQ[i,33] <- squarem$likelihood
    }

    test.qnem <- function(par){
      qnem <- QNEM(par, obs, modele.eye, M.step.eye, lower = rep(0,14), upper = rep(1,14), max.iter = Inf, trace.theta = FALSE)
      return(qnem)}
    d <- Sys.time()
    tr <- try(qnem <- test.qnem(par))
    f <- Sys.time()
    if(class(tr) != "try-error") {
      QNEM[i,1] <- qnem$iter
      QNEM[i,2:15] <- qnem$theta
      QNEM[i,16:29] <- par
      QNEM[i,30] <- qnem$forwards
      QNEM[i,31] <- qnem$backwards
      QNEM[i,32] <- f-d
      QNEM[i,33] <- qnem$neg.ll
    }
  }
  return(list(BW = BW , SQUAREM = SQ, qN = qN, QNEM = QNEM))
}


#observations
X.eye <- read.csv("/home/sidonie/Téléchargements/em-y35-fasttext.csv")
X.eye <- X.eye[,c(2,5,18)]
X.eye$READMODE[which(X.eye$READMODE < 0)] <- 0
X.eye$READMODE <- X.eye$READMODE +1
X.eye$PAIR <- paste0(X.eye$SUBJ, ":", X.eye$TEXT)
pairs <- unique(paste0(X.eye$SUBJ, ":", X.eye$TEXT))

#parameters
pi.eye <- p.stationnaire.eye(pi1 = 0.1, pi2 = 0.2)
par.eye <- c(  t12 <- 0.2,
               t13 <- 0.1,
               t21 <- 0.2,
               t23 <- 0.2,
               t31 <- 0.1,
               t32 <- 0.2,
               e12 <- 0.2,
               e13 <- 0.1,
               e21 <- 0.2,
               e23 <- 0.2,
               e31 <- 0.1,
               e32 <- 0.2,
               pi.eye[1],
               pi.eye[2])



#comparison test with 1000 different starting points
RNGkind("Mer") ; set.seed(28)
nb_it_eye <- test.eye(X.eye, par.eye, pairs)
saveRDS(nb_it_gd, "/home/sidonie/Bureau/results/HMM/res_gc.rds") #path to adapt to your environment
