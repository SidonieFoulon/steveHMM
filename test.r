require(steveHMM)
source("examples/umbrella_example.r")
em <- EM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, 200)
QNEM(par.parapluie, X.parapluie, modele.parapluie, M.step.parapluie, 200)

for(i in 1:50) {set.seed(i); QNEM(runif(2), X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 200)}

set.seed(7)
theta <- runif(2)
QNEM(theta, X.parapluie, modele.parapluie, M.step.parapluie, lower = c(0,0), upper = c(1,1), 200)
