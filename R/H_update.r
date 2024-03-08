# H + (s'y + y' H y) / (s'y)^2 * (ss') - (H y s' + s y' H) / s'y
H_update <- function(H, s, y) {
  sy <- sum(s * y)
  if(any(is.na(H))) { # As recommended in Nocedal Wright p 143
    H <- sy / sum(y**2) * diag(length(s))
  }
  Hy <- H %*% y
  yHy <- sum(y * Hy) 
  H + (sy + yHy)/sy**2 * tcrossprod(s,s) - (tcrossprod(Hy, s) + tcrossprod(s, Hy)) / sy
}
