
context("graph.complementer")

test_that("graph.complementer works", {

  library(igraph)

  g <- erdos.renyi.game(50, 3/50)
  g2 <- graph.complementer(g)
  g3 <- graph.complementer(g2)
  expect_that(graph.isomorphic(g, g3), is_true())
  
})


