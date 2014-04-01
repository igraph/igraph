
time_group("Fruchterman-Reingold layout")

time_that("FR layout is fast, connected", replications=10,
          init = { library(igraph); set.seed(42) },
          reinit = { g <- barabasi.game(400) },
          { layout.fruchterman.reingold(g, niter=500) })

time_that("FR layout is fast, unconnected", replications=10,
          init = { library(igraph); set.seed(42) },
          reinit = { g <- erdos.renyi.game(400, 400, type="gnm") },
          { layout.fruchterman.reingold(g, niter=500) })
