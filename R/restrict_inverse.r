
if(FALSE) {
# Vi = l'inverse d'une matrice symmétrique V 
# I = les indices des dimensions "bloquées"
restrict_inverse <- function(Vi, I) {
  V <- solve(Vi)
  V[I,] <- V[,I] <- 0
  R <- MASS::ginv(V)
  R[I,] <- R[,I] <- 0   # get rid of 1e-15's due to rounding errors
  R
}
}

# la formule est Vi - Vi %*% X %*% (t(X) %*% Vi %*% X)^{-1} %*% t(X)) %*% Vi
# où X engendre l'espace des dimensions bloquées
# dans notre cas X = [e_i : i \in I] et il y a des raccourcis...
restrict_inverse <- function(Vi, I) {
  n <- nrow(Vi)
  # ceci est (t(X) %*% Vi %*% X)^{-1}
  A <- solve( V[I,I] )
  # et ceci X %*% A %*% t(X)
  B <- matrix(0, nrow = n, ncol = n)
  B[I,I] <- A
  #
  R <- Vi - Vi %*% B %*% Vi
  R[I,] <- R[,I] <- 0   # get rid of 1e-15-s due to rounding errors
  R
}

