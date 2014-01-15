
time_group("SIR epidemics models on networks")

time_that("SIR is fast", replications=10,
          init = { library(igraph) },
          reinit = { g <- erdos.renyi.game(40, 40, type="gnm") },
          { sir(g, beta=5, gamma=1, no_sim=10) })
