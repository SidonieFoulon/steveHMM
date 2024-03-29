# uses salad::dual to get modele derivatives
modele_derivatives <- function(modele.fun, theta, obs) {
  theta1 <- dual(theta)
  vn <- varnames(theta1)

  modele <- modele.fun(theta1, obs)

  mod2 <- lapply(modele, function(x) if(is(x, "dual")) value(x) else x)
  mod2$dtrans <- d(modele$trans, vn)
  mod2$dp.emiss <- d(modele$p.emiss, vn)
  mod2$dpi <- d(modele$pi, vn)
  mod2
}
