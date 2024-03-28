# uses salad::dual to get modele derivatives
modele_derivatives <- function(modele.fun, theta, obs) {
  theta1 <- dual(theta)
  vn <- varnames(theta1)

  modele <- modele.fun(theta1, obs)

  list(trans = value(modele$trans),  p.emiss = value(modele$p.emiss),  pi = value(modele$pi),
      dtrans = d(modele$trans, vn), dp.emiss = d(modele$p.emiss, vn), dpi = d(modele$pi, vn) )
}
