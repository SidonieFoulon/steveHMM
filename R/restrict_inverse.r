#' @importFrom MASS ginv

if(TRUE) {
# Vi = l'inverse d'une matrice symmétrique V 
# I = les indices des dimensions "bloquées"
restrict_inverse <- function(Vi, I) {
  if(length(I) == 0) return(Vi)
  V <- MASS::ginv(Vi)
  V[I,] <- V[,I] <- 0
  R <- MASS::ginv(V)
  R[I,] <- R[,I] <- 0   # get rid of 1e-15's due to rounding errors
  R
}
} else {

# la formule est Vi - Vi %*% X %*% (t(X) %*% Vi %*% X)^{-1} %*% t(X)) %*% Vi
# où X engendre l'espace des dimensions bloquées
# dans notre cas X = [e_i : i \in I] et il y a des raccourcis...
# MAIS ça ne marche pas pour restreindre davantage une matrice déja projetée sur un sous ensemble de I
restrict_inverse <- function(Vi, I) {
  if(length(I) == 0) return(Vi)
  n <- nrow(Vi)
  # ceci est (t(X) %*% Vi %*% X)^{-1}
  A <- solve( Vi[I,I] )
  # et ceci X %*% A %*% t(X)
  B <- matrix(0, nrow = n, ncol = n)
  B[I,I] <- A
  #
  R <- Vi - Vi %*% B %*% Vi
  R[I,] <- R[,I] <- 0   # get rid of 1e-15's due to rounding errors
  R
}
}
