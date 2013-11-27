
time_group("local scan v1")

init   <- expression({library(igraph); set.seed(42) })
reinit <- expression({g <- random.graph.game(1000, p=.2)})

time_that("scan0 on random graphs",
          replications=10, init=init, reinit=reinit,
          { ls0 <- local.scan(g, k=0) })

time_that("scan1 on random graphs",
          replications=10, init=init, reinit=reinit,
          { ls1 <- local.scan(g) })
