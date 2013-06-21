
context("get.edgelist")

test_that("get.edgelist works", {
  library(igraph)
  g <- erdos.renyi.game(100, 3/100)
  e <- get.edgelist(g)
  g2 <- graph(t(e), n=vcount(g), dir=FALSE)
  expect_that(graph.isomorphic(g, g2), is_true())
})

