## on note : S in {Court : 1 ; Long : 2 ; Long stable : 3}
##           X in {Duree < 3min : 1 ; Duree >= 3min : 2}

modele.geyser.continu <- function(theta = c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39), obs, name.S = c("Court", "Long", "Long stable")) { #theta contient a, b + mu, sigma
  # parametre de la matrice de transition
  a <- theta[1]
  b <- theta[2]
  trans <- matrix(c(0,1-a,a,1,0,0,1-b,0,b), nrow = 3, byrow = TRUE)
  colnames(trans) <- rownames(trans) <- name.S

  # parametre de l'émission
  muc <- theta[3]
  mul <- theta[4]
  muls <- theta[5]
  sdc <- theta[6]
  sdl <- theta[7]
  sdls <- theta[8]

  #proba de l'observation par la loi normale
  p.emiss.c <- dnorm(obs, mean = muc, sd = sdc)
  p.emiss.l <- dnorm(obs, mean = mul, sd = sdl)
  p.emiss.ls <- dnorm(obs, mean = muls, sd = sdls)

  p.emiss <- rbind(p.emiss.c, p.emiss.l, p.emiss.ls)
  rownames(p.emiss) <- name.S

  # etat stationnaire
  pi_c <- theta[9]
  pi_l <- theta[10]
  pi_ls <- 1 - pi_c - pi_l
  pi <- c( c = pi_c, l = pi_l, ls = pi_ls )

  list(trans = trans, pi = pi, p.emiss = p.emiss)
}

M.step.geyser.continu <- function(obs, backward) {
  l <- ncol(backward$phi)

  #a = proba de changement d'état de "Court" à "Long stable" (->trans)
  D13 <- sum(backward$delta[1,3,-1])
  p.court <- sum(backward$phi[1,-l]) # = D12 + D13
  if(p.court > 0)
    a <- D13 / p.court
  else
    a <- 1 # pas d'état court

  #b = proba de passage d'état de "Long stable" à "Long stable" (->trans)
  D33 <- sum(backward$delta[3,3,-1])
  p.long.st <- sum(backward$phi[3,-l]) # = D31 + D33
  if(p.long.st > 0)
    b <- D33 / p.long.st
  else
    b <- 0 # pas d'état "long stable"

  #esperance de la loi normale pour les etats "court" et "Long", la proba pour l'état "Long stable" en découlera (->emiss)
  muc <- sum(backward$phi[1,] * obs) / sum(backward$phi[1,])

  mul <- sum(backward$phi[2,] * obs) / sum(backward$phi[2,])

  muls <- sum(backward$phi[3,] * obs) / sum(backward$phi[3,])

  #ecart-type de la loi normale pour les etats "court" et "Long", la proba pour l'état "Long stable" en découlera (->emiss)
  varc <- (sum(backward$phi[1,] * (obs**2)) / sum(backward$phi[1,])) - (muc**2)
  sdc <- sqrt(ifelse(varc > 0, varc, 0))

  varl <- (sum(backward$phi[2,] * (obs**2)) / sum(backward$phi[2,])) - (mul**2)
  sdl <- sqrt(ifelse(varl > 0, varl, 0))

  varls <- (sum(backward$phi[3,] * (obs**2)) / sum(backward$phi[3,])) - (muls**2)
  sdls <- sqrt(ifelse(varls > 0, varls, 0))

  #etats initiaux
  pi_c <- backward$phi[1,1]
  pi_l <- backward$phi[2,1]

  r <- c(a = a, b = b, muc = muc, mul = mul, muls = muls, sdc = sdc, sdl = sdl, sdls = sdls, pi_c = pi_c, pi_l = pi_l)
  r
}

# données geyser
X.geyser <- faithful$eruptions

# nos paramètres a et b et paramètres des lois normales d'initialisation :
par.geyser <- c(a = 0.31, b = 0.46, muc = 1.98, mul = 4.26, muls = 4.26, sdc = 0.28, sdl = 0.39, sdls = 0.39, pi_c = 0.44, pi_l = 0.30)



if(FALSE) {
library(steveHMM)

# qN
qn.cont <- quasi_newton(par.geyser, X.geyser, modele.geyser.continu, lower = c(.01,.01, rep(0, 6), .01, .01), upper = c(.99,.99, rep(Inf, 6), .99,.99))
qn.cont.trace <- capture_quasi_newton(par.geyser, X.geyser, modele.geyser.continu, lower = c(.01,.01, rep(0, 6), .01, .01), upper = c(.99,.99, rep(Inf, 6), .99,.99))

# EM
em.cont <- EM(par.geyser, X.geyser, modele.geyser.continu, M.step.geyser.continu, trace.theta = TRUE, criteria = "rel")

# SQUAREM
squarem.cont <- SQUAREM(par.geyser, X.geyser, modele.geyser.continu, M.step.geyser.continu, trace.theta = TRUE)

# QNEM
qnem.cont <- QNEM(par.geyser, X.geyser, modele.geyser.continu, M.step.geyser.continu,trace.theta = TRUE, verbose = TRUE)



#plots
par(mfrow = c(2,1))
plot( em$Theta[1,] , type = "o", xlim = c(0, max(ncol(em$Theta), ncol(squarem), ncol(qN))), ylim = range(em$Theta[1,], squarem[1,], qN[1,]), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter a")
abline(v = ncol(em$Theta), col = "#4080c0")
lines( squarem[1,], type = "o", col = "#c09140")
abline(v = ncol(squarem), col = "#c09140")
lines( qN[1,], type = "o", col = "#c05140")
abline(v = ncol(qN), col = "#c05140")
legend("topright", cex = 1, legend = c("Baum-Welch", "SQUAREM", "quasi-Newton"), col = c("#4080c0", "#c09140", "#c05140"), pch="o")

plot( em$Theta[2,] , type = "o", xlim = c(0, max(ncol(em$Theta), ncol(squarem), ncol(qN))), ylim = range(em$Theta[2,], squarem[2,], qN[2,]), col = "#4080c0", main = "Each values tested by the 3 algorithms for the parameter a")
abline(v = ncol(em$Theta), col = "#4080c0")
lines( squarem[2,], type = "o", col = "#c09140")
abline(v = ncol(squarem), col = "#c09140")
lines( qN[2,], type = "o", col = "#c05140")
abline(v = ncol(qN), col = "#c05140")
}





