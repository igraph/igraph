
time_group("Fruchterman-Reingold layout")

time_that("FR layout is fast", replications=10,
          init = { library(igraph); set.seed(42) },
          reinit = { g <- erdos.renyi.game(400, 400, type="gnm") },
          { layout.fruchterman.reingold(g, niter=500) })
