require(steveHMM)
source("eye_movement_trans.R")
X <- read.csv("em-y35-fasttext.csv")

X <- X[,c(2,5,18)]
X$READMODE[which(X$READMODE < 0)] <- 0
X$READMODE <- X$READMODE +1
X$PAIR <- paste0(X$SUBJ, ":", X$TEXT)

pairs <- unique(paste0(X$SUBJ, ":", X$TEXT))

OBS <- tapply(X$READMODE, X$PAIR, I)

par <- c(  t12 = 0.2, t13 = 0.1, t21 = 0.2, t23 = 0.2, t31 = 0.1, t32 = 0.2,
           e12 = 0.2, e13 = 0.1, e21 = 0.2, e23 = 0.2, e31 = 0.1, e32 = 0.2, 
           pi1 = 0.33, pi2 = 0.33 )


# contraintes lineaires
# tous sont positifs (C1 et U1), et on met les sommes 2 à 2 plus petites que 1 : ça implique
# que tous sont plus petit que 1
C1 <- -diag(14)
U1 <- rep(0,14)

#              t12 t13 t21 t23 t31 t32 e12 e13 e21 e23 e31 e32 pi1 pi2
C2 <- matrix(c(1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  
               0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
               0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  
               0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0, 
               0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
               0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0, 
               0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1), ncol = 14, byrow = TRUE)
U2 <- rep(1, 7)

C3 <- rbind(C1, C2)
U3 <- c(U1, U2)

# system.time( em <- EM(par, OBS, modele.eye, M.step.eye, max.iter=Inf, verbose = TRUE) )
# system.time(qnem <- QNEM2(par, OBS, modele.eye, M.step.eye, lin.coeff = C3, lin.upper = U3, max.iter = Inf, trace.theta = FALSE, verbose = TRUE))

random.par <- function() {
  ran <- function() {
    x <- runif(3)
    x[1:2] / sum(x)
  }

  r1 <- ran()
  r2 <- ran()
  r3 <- ran()
  r4 <- ran()
  r5 <- ran()
  r6 <- ran()
  r7 <- ran()

  par <- c(  t12 = r1[1], t13 = r1[2], 
             t21 = r2[1], t23 = r2[2], 
             t31 = r3[1], t32 = r3[2],
                          e12 = r4[1], 
                                        e13 = r5[1], 
             e21 = r6[1],               e23 = r5[2], 
             e31 = r6[2],  e32 = r4[2],
             pi1 = r7[1], pi2 = r7[2] )
  par
}

test.em <- function(seeds, data = OBS) {
  n <- length(seeds)
  D <- as.data.frame(rep(list(numeric(n)), 18))
  colnames(D) <- c("seed", "time", "iterations", "likelihood", "t12", "t13", "t21", "t23", 
                   "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2")
  for(i in seq_along(seeds)) {
    s <- seeds[i]
    set.seed(s)
    par <- random.par()
    time <- system.time( 
              tr <- try(
                em <- EM(par, data, modele.eye, M.step.eye, max.iter = 1000, verbose = FALSE, trace.theta = FALSE) 
              ) 
            )

    D$seed[i] <- s
    D$time[i] <- time[3]
    if(class(tr) != "try-error") {
      D$iterations[i] <- em$iter
      D$likelihood[i] <- -em$likelihood
      D[i, 5:18] <- em$theta
    }
    print(D[i,])
  }
  D
}


test.squarem <- function(seeds, data = OBS) {
  n <- length(seeds)
  D <- as.data.frame(rep(list(numeric(n)), 18))
  colnames(D) <- c("seed", "time", "iterations", "likelihood", "t12", "t13", "t21", "t23", 
                   "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2")
  for(i in seq_along(seeds)) {
    s <- seeds[i]
    set.seed(s)
    par <- random.par()
    time <- system.time( 
              tr <- try(
                squarem <- SQUAREM(par, data, modele.eye, M.step.eye, max.iter = 1000, trace.theta = FALSE) 
              ) 
            )

    D$seed[i] <- s
    D$time[i] <- time[3]
    if(class(tr) != "try-error") {
      D$iterations[i] <- squarem$iter
      D$likelihood[i] <- -squarem$likelihood
      D[i, 5:18] <- squarem$theta
    }
    print(D[i,])
  }
  D
}

test.qnem <- function(seeds, data = OBS) {
  n <- length(seeds)
  D <- as.data.frame(rep(list(numeric(n)), 20))
  colnames(D) <- c("seed", "time", "iterations", "fw", "bw", "likelihood", "t12", "t13", "t21", "t23", 
                   "t31", "t32", "e12", "e13", "e21", "e23", "e31", "e32", "pi1", "pi2")
  for(i in seq_along(seeds)) {
    s <- seeds[i]
    set.seed(s)
    par <- random.par()
    time <- system.time( 
              tr <- try(
                qnem <- QNEM2(par, data, modele.eye, M.step.eye, lin.coeff = C3, lin.upper = U3, max.iter = 1000, 
                              trace.theta = FALSE, verbose = FALSE)
              ) 
            )

    D$seed[i] <- s
    D$time[i] <- time[3]
    if(class(tr) != "try-error") {
      D$iterations[i] <- qnem$iter
      D$fw[i] <- qnem$forwards
      D$bw[i] <- qnem$backwards
      D$likelihood[i] <- qnem$neg.ll
     D[i, 7:20] <- qnem$theta
    }
    print(D[i,])
  }
  D
}

if(FALSE) {
Q <- test.qnem(1:6, OBS[1:50])
E <- test.em(1:6, OBS[1:50])

QNEM2(unlist(Q[1,7:20]), OBS[1:10],  modele.eye, M.step.eye, lin.coeff = C3, lin.upper = U3, max.iter = 1000, verbose = TRUE, trace.theta = FALSE)

set.seed(1)
par <- random.par()
QNEM2(par, OBS[1:10],  modele.eye, M.step.eye, lin.coeff = C3, lin.upper = U3, max.iter = 1000, verbose = TRUE, trace.theta = FALSE)
}


if(FALSE) {

for(p in pairs) {
  # obs <- X$READMODE[ X$PAIR == pairs[131] ]
  obs <- X$READMODE[ X$PAIR == p ]

  cat(p, "\n")
  if(length(unique(obs)) != 3) {
    cat("skipping! \n")
    next
  }

  # s <- system.time( EM(par, obs, modele.eye, M.step.eye, max.iter = Inf, trace.theta = FALSE, criteria = "rel") )
  s <- system.time( QNEM(par, obs, modele.eye, M.step.eye, lower = rep(0,14), upper = rep(1,14), max.iter = Inf, trace.theta = FALSE) )
  print(s)
}

}

