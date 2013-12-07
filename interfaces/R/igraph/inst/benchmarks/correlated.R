
time_group("correlated E-R graphs, v1")

time_that("correlated.game is fast", replications=10,
          init={ library(Matrix); library(igraph) },
          { correlated.pair.game(100, corr=.8, p=5/100) })


