#' @importFrom salad dual
#' @importFrom salad varnames
#' @importFrom salad value
#' @importFrom salad d
#' @importFrom methods is

modele_derivatives <- function(modele.fun, theta, obs) {
  theta1 <- dual(theta)
  vn <- varnames(theta1)

  modele <- modele.fun(theta1, obs)

  mod2 <- lapply(modele, function(x) if(is(x, "dual")) value(x) else x)
 
  mod2$dtrans <- d(modele$trans, vn)

  if(is.list(modele$p.emiss)) {
    mod2$p.emiss <- lapply(modele$p.emiss, value)
    mod2$dp.emiss <- lapply(modele$p.emiss, d, varnames = vn)
  } else {
    mod2$dp.emiss <- d(modele$p.emiss, vn)
  }

  mod2$dpi <- d(modele$pi, vn)
  mod2
}
