
context("largest.cliques")

test_that("largest.cliques works", {

  library(igraph)

  g <- erdos.renyi.game(50,20/50)
  lc <- largest.cliques(g)

  ## TODO: this only checks that these are cliques
  expect_that(unique(sapply(lc, function(x)
                            graph.density(induced.subgraph(g, x)))),
              equals(1))

})
