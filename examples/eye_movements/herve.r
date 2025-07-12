require(steveHMM)
source("eye_movement_trans.R")
X <- read.csv("em-y35-fasttext.csv")

X <- X[,c(2,5,18)]
X$READMODE[which(X$READMODE < 0)] <- 0
X$READMODE <- X$READMODE +1
X$PAIR <- paste0(X$SUBJ, ":", X$TEXT)

pairs <- unique(paste0(X$SUBJ, ":", X$TEXT))

L <- tapply(X$READMODE, X$PAIR, I)

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
C <- rbind(C1, C2)
U <- c(U1, U2)

system.time(qnem <- QNEM2(par, L, modele.eye, M.step.eye, lin.coeff = C, lin.upper = U, max.iter = Inf, trace.theta = FALSE, verbose = TRUE))


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


system.time( em <- EM(par, L, modele.eye, M.step.eye, max.iter=Inf, verbose = TRUE) )
system.time(qnem <- QNEM(par, L, modele.eye, M.step.eye, lower = rep(0,14), upper = rep(1,14), max.iter = Inf, trace.theta = FALSE, verbose = TRUE))


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
C <- rbind(C1, C2)
U <- c(U1, U2)

system.time(qnem <- QNEM2(par, L, modele.eye, M.step.eye, lin.coeff = C, lin.upper = U, max.iter = Inf, trace.theta = FALSE, verbose = TRUE))


stop()
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

# p <- "5:associations_humanitaires-a1"

obs <- X$READMODE[ X$PAIR == p ]
QNEM(par, obs, modele.eye, M.step.eye, lower = rep(0,14), upper = rep(1,14), max.iter = Inf, trace.theta = FALSE, verbose = TRUE)


