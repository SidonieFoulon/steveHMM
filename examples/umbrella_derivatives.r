## on note : S in {Dry : 1 ; Rainy : 2}
##           X in {No umbrella : 1 ; Umbrella : 2}

# le modèle avec les dérivées...
modele.parapluie2 <- function(theta, obs, name.S = c("Dry", "Rainy")) {
  a <- theta[1]
  b <- theta[2]

  # matrice de transition
  trans <- matrix( c(1-a, a, a, 1-a), byrow = TRUE, nrow = 2, ncol = 2)
  trans.da <- matrix( c(-1, 1, 1, -1), byrow = TRUE, nrow = 2, ncol = 2)
  trans.db <- matrix(0, nrow = 2, ncol = 2)
  colnames(trans) <- rownames(trans) <- name.S

  # etat stationnaire, ici ne dépend pas de a ou b
  stat <- c(0.5, 0.5)
  stat.da <- c(0, 0)
  stat.db <- c(0, 0)

  # probas d'emission
  p.emiss <- rbind(ifelse(obs == 1, 1-b, b), ifelse(obs == 1, b, 1-b))
  p.emiss.da <- matrix(0, nrow = 2, ncol = length(obs))
  p.emiss.db <- rbind(ifelse(obs == 1, -1, 1), ifelse(obs == 1, 1, -1))
  rownames(p.emiss) <- name.S

  list(trans = trans, pi = stat, p.emiss = p.emiss, 
       d.trans = list(trans.da, trans.db), d.pi = list(stat.da, stat.db),
       d.p.emiss = list(p.emiss.da, p.emiss.db) )
}
