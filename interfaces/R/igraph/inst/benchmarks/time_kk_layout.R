
time_group("Kamada-Kawai layout")

time_that("KK layout is fast, connected", replications=10,
          init = { library(igraph); set.seed(42) },
          reinit = { g <- barabasi.game(400) },
          { layout.kamada.kawai(g, maxiter=500) })

time_that("KK layout is fast, unconnected", replications=10,
          init = { library(igraph); set.seed(42) },
          reinit = { g <- erdos.renyi.game(400, 400, type="gnm") },
          { layout.kamada.kawai(g, maxiter=500) })

time_that("KK layout is fast for large graphs", replications=10,
          init = { library(igraph); set.seed(42) },
          reinit = { g <- barabasi.game(3000) },
          { layout.kamada.kawai(g, maxiter=500) })
