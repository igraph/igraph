
context("largest.independent.vertex.sets")

test_that("largest.independent.vertex.sets works", {

  library(igraph)

  g <- erdos.renyi.game(50, 0.8)
  livs <- largest.independent.vertex.sets(g)
  expect_that(unique(sapply(livs, length)),
              equals(independence.number(g)))

  ec <- sapply(seq_along(livs), function(x)
               ecount(induced.subgraph(g, livs[[x]])))
  expect_that(unique(ec), equals(0))

  ## TODO: check that they are largest
})
